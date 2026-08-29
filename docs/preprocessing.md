# Trimming and filtering

The `trim` module is the fastp-lite half of preprocessing: adapter and quality
trimming, length and content filtering, and subsampling. Every entry point
streams, and output is gzipped automatically when the target name ends in `.gz`.

```jennifer
import "@jennifer/ngs/trim.j" as trim;

def opts as trim.Options init trim.options();
$opts.qualityThreshold = 20;
$opts.minLength = 36;

def st as trim.Stats init trim.run("in.fastq.gz", "out.fastq.gz", $opts);
io.printf("%d of %d reads kept\n", $st.readsOut, $st.readsIn);
```

## The order of operations

`trim.run` applies the same steps fastp does, in the same order. Each is also
exposed on its own, so a pipeline can apply only the parts it wants.

1. **fixed clipping** - `clipFixed(rec, front, tail)`
2. **adapter removal** - `clipAdapters(rec, adapters, minOverlap)`
3. **quality trimming** - `clipQuality(rec, threshold, window, offset)`
4. **length capping** - `truncate(rec, maxLength)`

Then filtering, on the trimmed record: `keeps(rec, opts)`.

`trim.trimRead(rec, opts)` runs all four steps and returns the result, for when
you want the transformation without the streaming.

## Options

The zero value is a no-op: it trims nothing and keeps everything.
`trim.options()` returns useful defaults - adapter clipping on, sliding-window
quality trimming at Q20, reads under 25 bases discarded.

| field | default from `options()` | meaning |
| --- | --- | --- |
| `trimFront` / `trimTail` | `0` | bases removed from each end |
| `adapters` | built-in set | adapters to clip; empty disables clipping |
| `minOverlap` | `6` | shortest partial adapter accepted at the 3' end |
| `qualityThreshold` | `20` | window mean quality below which the read is cut; `0` disables |
| `windowSize` | `4` | sliding window width |
| `maxLength` | `0` | truncate longer reads; `0` is no cap |
| `minLength` | `25` | discard reads shorter than this after trimming |
| `maxAmbiguous` | `-1` | discard reads with more N bases than this; `-1` disables |
| `minMeanQuality` | `0` | discard reads below this mean quality; `0` disables |
| `offset` | `PHRED33` | quality encoding |

## Adapter clipping

Two cases, and both matter.

A **full match** anywhere in the read is found with `strings.indexOf` - a native
substring scan. The read is cut at the adapter's start.

A **partial adapter running off the 3' end** is the common case in real data:
the fragment was shorter than the read, so the read ends in the first few bases
of the adapter and never reaches the rest of it. Those are found by testing
adapter prefixes against the read's suffix, longest first, down to `minOverlap`
bases.

```jennifer
trim.adapterStart($seq, qc.defaultAdapters(), 6);   # cut position, or -1
```

When several adapters match, the earliest cut wins.

`minOverlap` is a real trade-off. Six bases of random sequence match a given
adapter prefix about once in 4 096 reads, so a low value trims a few clean reads
by a handful of bases. Raising it misses genuine short overhangs. Six is the
usual compromise; real trimmers behave the same way.

## Quality trimming

The Trimmomatic and fastp sliding-window rule: scan windows from the 5' end and
cut at the first window whose mean quality falls below the threshold. This
handles the normal degradation pattern, where quality falls off toward the 3'
end.

```jennifer
trim.qualityKeepLength("IIII####", 20, 2, 33);   # 4 - keeps the good prefix
```

The running window sum is maintained incrementally, so this is one pass over the
quality string, not one per window. It is still per-base work - see
[Performance](performance.md).

## Filtering

`keeps` applies the three content filters to an already-trimmed record:
minimum length, maximum ambiguous bases, minimum mean quality. `trim.run`
applies them for you and records *why* each read was dropped.

```jennifer
def st as trim.Stats init trim.run($in, $out, $opts);

$st.readsIn; $st.readsOut;
$st.basesIn; $st.basesOut;
$st.adapterTrimmed;    # records shortened by adapter clipping
$st.qualityTrimmed;    # records shortened by quality trimming
$st.tooShort;
$st.tooAmbiguous;
$st.lowQuality;
```

## Subsampling

Two ways, both streaming.

**Deterministic**, and the one to reach for - reproducible with no seed to
record, and the natural way to cut a large run down to a test-sized sample:

```jennifer
trim.subsampleEvery("in.fastq.gz", "out.fastq.gz", 10, 0);       # every 10th read
trim.subsampleEvery("in.fastq.gz", "out.fastq.gz", 1, 100000);   # the first 100k
```

**Random**, when you need an unbiased sample rather than a regular one. Seed
`math.randSeed` first to make it reproducible:

```jennifer
use math;
math.randSeed(42);
trim.subsampleFraction("in.fastq.gz", "out.fastq.gz", 0.1);      # ~10% of reads
```

Both return a `trim.Stats` with the read and base counts filled in.
