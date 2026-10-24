# Trimming and filtering

`trim` does adapter and quality trimming, length and content filtering, and
subsampling, in the manner of fastp. Every entry point streams, and output is
gzipped when the target name ends in `.gz`.

```jennifer
import "@deckyard/ngs/trim.j" as trim;

def opts as trim.Options init trim.options();
$opts.qualityThreshold = 20;
$opts.minLength = 36;

def st as trim.Stats init trim.run("in.fastq.gz", "out.fastq.gz", $opts);
io.printf("%d of %d reads kept\n", $st.readsOut, $st.readsIn);
```

`trim.run`, `trim.subsampleEvery` and `trim.subsampleFraction` raise
`ngs.trim` when the input and output are the same file (symlinks resolved):
the output is opened for writing before the input is read and would be
truncated. They also raise `ngs.trim` when the compressor exits non-zero. The
input reader is closed with `close`, not `closeStrict`, so a truncated `.gz`
input ends the run early without an error.

## The order of operations

`trim.run` applies fastp's steps in fastp's order. Each is exported, so a
program can apply only some of them:

1. fixed clipping - `trim.clipFixed(rec, front, tail)`
2. adapter removal - `trim.clipAdapters(rec, adapters, minOverlap)`
3. quality trimming - `trim.clipQuality(rec, threshold, window, offset)`
4. length capping - `trim.truncate(rec, maxLength)`

Filtering, `trim.keeps(rec, opts)`, then runs on the trimmed record.
`trim.trimRead(rec, opts)` applies the four steps to one record without
streaming.

## Options

The zero value trims nothing and applies no filter except `maxAmbiguous` (see
below the table). `trim.options()` turns on
adapter clipping, sliding-window quality trimming at Q20, and discards reads
under 25 bases.

| field | from `options()` | meaning |
| --- | --- | --- |
| `trimFront` / `trimTail` | `0` | bases removed from the 5' / 3' end |
| `adapters` | `qc.defaultAdapters()` | adapters to clip; empty disables clipping |
| `minOverlap` | `6` | shortest partial adapter accepted at the 3' end; `0` means `trim.DEFAULT_OVERLAP` (6) |
| `qualityThreshold` | `20` | window mean quality below which the read is cut; `0` disables |
| `windowSize` | `4` | sliding window width; `0` means `trim.DEFAULT_WINDOW` (4) |
| `maxLength` | `0` | truncate longer reads; `0` is no cap |
| `minLength` | `25` | discard reads shorter than this after trimming |
| `maxAmbiguous` | `-1` | discard reads with more N bases than this; `-1` disables |
| `minMeanQuality` | `0` | discard reads below this mean quality; `0` disables |
| `offset` | `fastq.PHRED33` | quality encoding; `0` means Phred+33 |

The zero value's `maxAmbiguous` is `0`, which discards any read containing an
N; set it to `-1` when building `Options` by hand.

## Adapter clipping

A full adapter anywhere in the read is found with `strings.indexOf`, a native
substring scan, and the read is cut at its start.

A partial adapter running off the 3' end is the usual case: the fragment was
shorter than the read, which ends in the first bases of the adapter. These are
found by testing adapter prefixes against the read's suffix, longest first,
down to `minOverlap` bases.

```jennifer
trim.adapterStart($seq, qc.defaultAdapters(), 6);   # 0-based cut position, or -1
```

With several adapters the earliest cut wins.

Six random bases match a given adapter prefix about once in 4 096 reads, so a
low `minOverlap` shortens some clean reads by a few bases, and a high one
misses short overhangs. Six is the usual compromise.

## Quality trimming

The Trimmomatic / fastp sliding-window rule: scan windows from the 5' end and
cut at the start of the first window whose mean quality falls below the
threshold.

```jennifer
trim.qualityKeepLength("IIII####", 20, 2, 33);   # 4 - keeps the good prefix
```

A window wider than the read shrinks to the read. The window sum is updated
incrementally, so this is one pass over the quality string, still one
interpreter step per base ([Performance](performance.md)).

## Filtering

`trim.keeps` applies minimum length, maximum ambiguous bases and minimum mean
quality to a trimmed record. `trim.run` applies the same filters and records
why each read was dropped:

```jennifer
def st as trim.Stats init trim.run($in, $out, $opts);

$st.readsIn; $st.readsOut;
$st.basesIn; $st.basesOut;
$st.adapterTrimmed;    # records shortened by adapter clipping
$st.qualityTrimmed;    # records shortened by quality trimming
$st.tooShort;          # dropped by minLength
$st.tooAmbiguous;      # dropped by maxAmbiguous
$st.lowQuality;        # dropped by minMeanQuality
```

## Subsampling

Deterministic, reproducible without a seed:

```jennifer
trim.subsampleEvery("in.fastq.gz", "out.fastq.gz", 10, 0);       # every 10th read
trim.subsampleEvery("in.fastq.gz", "out.fastq.gz", 1, 100000);   # the first 100k
```

An `everyNth` of 1 or less keeps every read; a `maxReads` of 0 means no limit.

Random, for an unbiased rather than regular sample. Seed `math.randSeed` first
to make it reproducible:

```jennifer
use math;
math.randSeed(42);
trim.subsampleFraction("in.fastq.gz", "out.fastq.gz", 0.1);      # ~10% of reads
```

The fraction must be in (0, 1], otherwise `ngs.trim` is raised. Both return a
`trim.Stats` with only the read and base counts filled in.
