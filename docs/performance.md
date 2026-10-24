# Performance

The tree-walking interpreter costs roughly a microsecond per statement. Line
I/O is cheap; loops that touch every base are what make Jennifer code slow.

## The measurements

50 000 reads of 150 bp (7.5 M bases, a 15 MB uncompressed FASTQ), one core:

| operation | time | throughput |
| --- | --- | --- |
| read lines, convert quality to `bytes`, count | 0.28 s | ~54 MB/s |
| naive per-base loop over the same bases | 14.1 s | ~0.5 M bases/s |
| `qc.analyze`, default options | 3.1 s | ~16 100 reads/s |
| `qc.analyze` with `perBase: true` | 27.5 s | ~1 800 reads/s |

`fs.readLine` and `convert.bytesFromString` run in Go. The 50x gap between the
first two rows is the cost of the interpreter loop.

## Count in Go

`core.countIn(haystack, needle)` is:

```jennifer
(len($haystack) - len(strings.replace($haystack, $needle, ""))) // len($needle)
```

`strings.replace` runs in Go, so this is one native scan. [`qc`](qc.md) buffers
a batch of records and then counts:

- GC and N content - `core.gcCount` and `core.nCount`, four and two scans per
  batch
- adapters - one scan per probe per batch
- the quality histogram - one scan per Phred score, stopping once every
  character in the batch is accounted for: about forty scans for typical
  Illumina data instead of the whole printable-ASCII range

The counts are exact.

## What is still per-base

| operation | why it cannot batch |
| --- | --- |
| `qc` with `perBase: true` | position *i* of every read is not a character count |
| `trim.qualityKeepLength` | a sliding window is positional |
| `fastq.qualities` / `meanQuality` / `minQuality` | each score is decoded individually |

150 steps per read is negligible for a few reads and significant across 150
million.

## Budgeting a run

At ~16 100 reads/s, a full QC pass over a 30x human WGS run of about 150 M reads
takes around three hours. Subsampling takes seconds and gives stable
statistics:

```jennifer
def opts as qc.Options init qc.options();
$opts.everyNth = 20;         # spread across the whole file
$opts.maxReads = 200000;     # and stop at 200k
```

Trimming has no such shortcut, since every read has to be written. For a large
run, have [`pipeline`](pipelines.md) call `fastp` and parse its report.

## Memory is bounded by the batch

Streaming keeps the file out of memory, but the QC batch is held in memory:

| input | batch bound | peak RSS |
| --- | --- | --- |
| 50 000 x 150 bp | 2000 records (300 KB) | 34 MB |
| 3 000 x 30 kb, records only | 2000 records (120 MB) | 657 MB |
| 3 000 x 30 kb, records + bases | 4 MB of bases | 74 MB |

The third row is `qc.Options.maxBatchBases` at its default; results are the
same either way. When buffering records for Go, bound the buffer in bytes:
read lengths span three orders of magnitude.

## Struct copies in hot loops

Assignment and argument passing copy. In a module, a parameter the body never
writes is borrowed, not copied, so passing a large interval list to
`intervals.intersect` is free; a struct built up per record is not.

An early `qc.analyze` accumulated into the `Report` struct, deep-copying a
94-element histogram and two maps per record. Moving the hot counters into
locals and building the `Report` at the end took the default path from 6.5 s to
4.0 s at the time. In a per-record loop, accumulate into scalars and locals;
batch-level state can live in a struct, copied 2000x less often.

## Building lists

```jennifer
$xs[] = item;                      # amortised O(n) for n appends
$xs = lists.push($xs, item);       # O(n^2) - copies the whole list every pass
```

Two related traps the deck hit:

- `$m[k][] = v` is a parse error, and the read-modify-write workaround is
  quadratic. `intervals.byChromosome` sorts and groups runs in one pass
  instead.
- `binary.concat` in a loop is O(n^2) too. The deck accumulates a
  `list of string` and joins once.

## Parallelism

`pipeline.runParallel` fans external processes across cores; that is where
throughput comes from, since the work runs in optimised native code.

`spawn` runs Jennifer concurrently. It deep-copies its enclosing scope at
launch, so there are no shared accumulators. It suits per-sample work that
returns a small summary; it does not make a per-base loop fast.

## Profiler findings

`jennifer profile` reports hits and wall-clock time per source position;
`--allocs` reports every deep copy. Both were run against 200 000-record
fixtures while tuning the deck.

### Binding a list element copies the struct

```jennifer
for (def i in 1..len($ivs)) {
    def current as Interval init $ivs[$i];   # deep-copies the struct, every pass
    if ($current.start < $limit) { ... }
}
```

One such line in `intervals.merge` caused 400 000 copies. Reading fields off
the indexed element means the same and allocates nothing:

```jennifer
for (def i in 1..len($ivs)) {
    if ($ivs[$i].start < $limit) { ... }
}
```

`for (def iv in $xs)` over a list of structs also copies per element. It
matters on per-record paths, not over small lists or lists of ints.

### One exception frame per record

Measured: `convert.toInt` 1.4 µs, wrapped in a function 2.5 µs, wrapped in a
function with a `try` 3.5 µs, so a `try` per column costs about 2 µs. SAM has
five numeric columns. Wrapping the whole record parse in one `try` and
re-raising through `core.reraise` costs one frame, and the message names the
record rather than the column.

### Helper calls in the innermost loop

`core.field($f, 3, "")` at two calls per record was measurably slower than an
inline `if (len($f) > 3)`. The format parsers inline their column reads, with a
comment saying why.

### Check before sorting

`intervals.sort` builds one zero-padded key per element for `lists.sortBy`. BED
and GFF files usually arrive sorted, and every set operation sorts its inputs,
so `sort` first scans for order, comparing fields with no allocation, and stops
at the first pair out of order. On a 200 000-feature BED, read-merge-measure
takes 16.9 s shuffled and 12.5 s sorted, and the sorted run allocates 480 MB
less.

## Where to spend effort

1. Read less: subsample, or filter early in the stream.
2. Count in Go: `core.countIn`, `strings.indexOf`, `strings.split`.
3. Do not copy in the hot loop: index list elements, keep per-record
   accumulators in scalars.
4. Profile first and measure after. Three changes that looked obviously good
   did not survive measurement: a cheaper sort key bought 9%, building the key with
   `strings.join` was slower than concatenation, and a single-character
   pre-filter in the adapter matcher, skipping a slice-and-compare three times in
   four, was 18% slower because its own slice cost as much as the one it
   avoided. All three were reverted. Doing less work and running fewer
   statements are different claims, and only the second is reliably faster
   here.
5. Use the native tool when the work is per-base over a whole run;
   [`pipeline`](pipelines.md) drives it.
