# Read quality control

`qc.analyze` makes a single streaming pass over a FASTQ file and returns
everything it measured: read and base counts, the length distribution, GC and N
content, the full quality-score histogram, adapter occurrence counts, and
optionally per-position mean quality.

```jennifer
use io;
import "@jennifer/ngs/qc.j" as qc;

def rep as qc.Report init qc.analyze("sample.fastq.gz", qc.options());
io.printf("%s", qc.summary($rep));
```

```
reads seen      4000
reads sampled   4000
bases           600000
read length     150-150 (mean 150.0)
GC content      50.09 %
N bases         0
mean quality    32.13
bases >= Q30    69.57 %
adapter AGATCGGAAGAGC 800
adapter CTGTCTCTTATACACATCT 0
adapter TGGAATTCTCGG 0
adapter GGGGGGGGGG 1
```

## Options

`qc.Options` has a valid zero value - every read, no per-position quality,
Phred+33 - so you can enable one thing at a time. `qc.options()` returns the
same defaults explicitly.

| field | default | meaning |
| --- | --- | --- |
| `maxReads` | `0` | stop after sampling this many reads; `0` is no limit |
| `everyNth` | `1` | sample one read in every N |
| `perBase` | `false` | also collect per-position quality |
| `offset` | `PHRED33` | quality encoding |
| `batchSize` | `2000` | records buffered per counting batch |
| `maxBatchBases` | `4000000` | bases buffered per batch, whichever limit fires first |
| `adapters` | built-in set | sequences to count; empty means the built-ins |

```jennifer
def opts as qc.Options init qc.options();
$opts.maxReads = 200000;
$opts.perBase = true;
def rep as qc.Report init qc.analyze("sample.fastq.gz", $opts);
```

## Derived statistics

The `Report` holds raw counts; these turn them into the numbers you report.

```jennifer
qc.meanQuality($rep);    # mean Phred across every sampled base
qc.q30Percent($rep);     # percentage of bases at Q30 or better
qc.gcPercent($rep);      # percentage G+C
qc.meanLength($rep);     # mean sampled read length
qc.perBaseMean($rep);    # list of float, mean quality per position
qc.summary($rep);        # the text block above
```

Everything else is on the struct: `reads`, `sampled`, `bases`, `minLength`,
`maxLength`, `gcBases`, `nBases`, `qualHist` (indexed by Phred score),
`lengthHist`, `adapterHits`.

## Adapter detection

Four probes are counted by default:

| constant | sequence | what it is |
| --- | --- | --- |
| `qc.ADAPTER_TRUSEQ` | `AGATCGGAAGAGC` | Illumina TruSeq / NEBNext universal |
| `qc.ADAPTER_NEXTERA` | `CTGTCTCTTATACACATCT` | Nextera / transposase |
| `qc.ADAPTER_SMALL_RNA` | `TGGAATTCTCGG` | Illumina small-RNA 3' |
| `qc.ARTEFACT_POLY_G` | `GGGGGGGGGG` | two-colour chemistry poly-G artefact |

`qc.defaultAdapters()` returns that list; set `Options.adapters` to override it.

These are **occurrence counts**, not read counts - a read containing the adapter
twice contributes two. In practice the two are nearly identical, since an
adapter appears once per read.

Sequences are joined with newlines before scanning, so an adapter can never be
matched across a record boundary. That is a real hazard of the batching approach
and the test suite pins it.

## Why it is not slow

A tree-walking interpreter costs roughly a microsecond per step, so touching
every base in Jennifer is what makes naive QC slow: a per-base loop over 7.5 M
bases takes about 14 seconds.

`qc.analyze` instead buffers a batch of records and counts characters with
`core.countIn`, which is `strings.replace` inside Go - one native scan per
character class per batch, rather than one interpreter step per base. The
quality histogram is built the same way: one pass per Phred score, stopping as
soon as every character in the batch is accounted for, so a typical Illumina
batch costs about forty native passes rather than the full printable range.

The results are identical, not approximate; they are checked against an
independent reference implementation in the test suite.

A batch is bounded by **both** a record count and a base count. The record
count alone is only a memory bound when reads are short: 2000 nanopore reads of
30 kb are 60 MB of sequence plus 60 MB of quality before the joins, which
measured 657 MB resident. `maxBatchBases` caps that at 74 MB for the same file
without changing any result, and never fires for short reads, where the record
limit is still the one that applies.

See [Performance](performance.md) for the measured numbers.

## Per-position quality

Per-position mean quality is the one statistic that cannot be expressed as a
character count - you have to look at position *i* of every read - so it costs
one interpreter step per base and is opt-in.

```jennifer
def opts as qc.Options init qc.options();
$opts.perBase = true;
$opts.maxReads = 50000;          # pair it with subsampling

def rep as qc.Report init qc.analyze("sample.fastq.gz", $opts);
def means as list of float init qc.perBaseMean($rep);
```

Ragged read lengths are handled: `perBaseCount` records how many reads reached
each position, so the mean at position 200 is over only the reads that were that
long.

## Subsampling

For a large run, sample. This is what FastQC does too, and the statistics are
stable long before the file is.

```jennifer
$opts.maxReads = 200000;   # the first 200k reads
$opts.everyNth = 20;       # every 20th read, spread across the whole file
```

`everyNth` is the better sample when a run's quality drifts over its length,
since it spreads across the file rather than reading only the front. The two
compose: `everyNth = 20` with `maxReads = 200000` samples every 20th read until
it has 200 000 of them.

`Report.reads` counts the records read, `Report.sampled` the records measured.

## Wrong encoding

If a quality character decodes to a score outside 0-93, `analyze` raises
`ngs.qc` naming the likely cause - an offset mismatch. Legacy Illumina 1.3-1.7
data needs `$opts.offset = fastq.PHRED64;`.
