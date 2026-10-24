# Read quality control

`qc.analyze` makes one streaming pass over a FASTQ file and returns read and
base counts, the length distribution, GC and N content, the quality-score
histogram, adapter occurrence counts and, on request, per-position mean
quality.

```jennifer
use io;
import "@deckyard/ngs/qc.j" as qc;

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

The zero value of `qc.Options` is valid: every read, no per-position quality,
Phred+33. `qc.options()` returns the same settings explicitly.

| field | default | meaning |
| --- | --- | --- |
| `maxReads` | `0` | stop after sampling this many reads; `0` is no limit |
| `everyNth` | `1` | sample one read in every N; `0` or `1` is every read |
| `perBase` | `false` | also collect per-position quality |
| `offset` | `fastq.PHRED33` | quality encoding; `0` means Phred+33 |
| `batchSize` | `2000` (`qc.DEFAULT_BATCH`) | records buffered per counting batch |
| `maxBatchBases` | `4000000` (`qc.DEFAULT_BATCH_BASES`) | bases buffered per batch; whichever limit is reached first counts the batch |
| `adapters` | empty | sequences to count; empty means `qc.defaultAdapters()` |

```jennifer
def opts as qc.Options init qc.options();
$opts.maxReads = 200000;
$opts.perBase = true;
def rep as qc.Report init qc.analyze("sample.fastq.gz", $opts);
```

`analyze` raises `ngs.fastq` on a malformed record and `ngs.source` on a
missing file.

## Derived statistics

```jennifer
qc.meanQuality($rep);    # mean Phred over every sampled base
qc.q30Percent($rep);     # percentage of bases at Q30 or better
qc.gcPercent($rep);      # percentage G+C
qc.meanLength($rep);     # mean sampled read length
qc.perBaseMean($rep);    # list of float, mean quality per position
qc.summary($rep);        # the text block above
```

The four scalar functions return `0.0` when nothing was sampled;
`perBaseMean` returns an empty list unless `perBase` was set. The raw counts
are fields of `qc.Report`: `reads`, `sampled`, `bases`, `minLength`,
`maxLength`, `gcBases`, `nBases`, `qualHist` (indexed by Phred score, 0 to 93),
`lengthHist`, `adapterHits`, and `perBaseSum` / `perBaseCount` when `perBase` is
set.

## Adapter detection

Four probes are counted by default:

| constant | sequence | source |
| --- | --- | --- |
| `qc.ADAPTER_TRUSEQ` | `AGATCGGAAGAGC` | Illumina TruSeq / NEBNext universal |
| `qc.ADAPTER_NEXTERA` | `CTGTCTCTTATACACATCT` | Nextera / transposase |
| `qc.ADAPTER_SMALL_RNA` | `TGGAATTCTCGG` | Illumina small-RNA 3' |
| `qc.ARTEFACT_POLY_G` | `GGGGGGGGGG` | poly-G artefact of two-colour chemistry (NovaSeq, NextSeq) |

Counts are occurrences, not reads: a read holding the adapter twice adds two.
Sequences in a batch are joined with newlines before scanning, so a match never
spans two records; a test covers this.

## How it counts

`qc.analyze` buffers a batch of records and counts characters with
`core.countIn`, a `strings.replace` inside Go: one native scan per character
class per batch rather than one interpreter step per base. The quality
histogram takes one scan per Phred score and stops once every character in the
batch is accounted for, about forty scans for typical Illumina data. Results are
exact; the tests check them against hand-computed values and show that moving a
batch boundary does not change them.

A batch is bounded by a record count and a base count. The record count alone
bounds memory only for short reads: 2000 nanopore reads of 30 kb measured
657 MB resident, and `maxBatchBases` brings the same file to 74 MB with the same
results. Timings and memory figures are in [Performance](performance.md).

## Per-position quality

Mean quality at position *i* cannot be expressed as a character count, so
`perBase` costs one interpreter step per base. Pair it with subsampling:

```jennifer
def opts as qc.Options init qc.options();
$opts.perBase = true;
$opts.maxReads = 50000;

def rep as qc.Report init qc.analyze("sample.fastq.gz", $opts);
def means as list of float init qc.perBaseMean($rep);
```

Reads of different lengths are handled: `perBaseCount` records how many reads
reached each position, and the mean at a position is over those reads only.

## Subsampling

The statistics settle long before the end of a large run; FastQC samples too.

```jennifer
$opts.maxReads = 200000;   # stop after 200k sampled reads
$opts.everyNth = 20;       # every 20th read, spread across the file
```

`everyNth` samples the whole file, so it catches quality that drifts over a
run; `maxReads` alone reads only the front. Together, `everyNth = 20` and
`maxReads = 200000` take every 20th read until 200 000 are sampled.
`Report.reads` counts records read, `Report.sampled` records measured.

## Wrong encoding

A quality character that decodes outside Phred 0-93 raises `ngs.qc`, naming an
offset mismatch as the likely cause. Legacy Illumina 1.3-1.7 data needs
`$opts.offset = fastq.PHRED64;`.
