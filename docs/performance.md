# Performance

One fact explains every performance decision in this deck:

> A tree-walking interpreter costs roughly a **microsecond per step**. Touching
> every base in Jennifer is what makes bioinformatics code slow.

Line I/O is not the problem. Per-base loops are.

## The measurements

All figures below are from 50 000 × 150 bp reads — 7.5 M bases, a 15 MB
uncompressed FASTQ — on one core.

| operation | time | throughput |
| --- | --- | --- |
| read lines, convert quality to `bytes`, count | 0.28 s | ~54 MB/s |
| naive per-base loop over the same bases | 14.1 s | ~0.5 M bases/s |
| `qc.analyze`, default options | 3.1 s | ~16 100 reads/s |
| `qc.analyze` with `perBase: true` | 27.5 s | ~1 800 reads/s |

`fs.readLine` and `convert.bytesFromString` are Go, so getting data *in* is
cheap. The 50× gap between rows one and two is the cost of the interpreter loop
itself.

## The technique: count in Go, not in Jennifer

`core.countIn(haystack, needle)` is:

```jennifer
(len($haystack) - len(strings.replace($haystack, $needle, ""))) // len($needle)
```

`strings.replace` runs inside Go, so this is one native scan rather than a
per-character loop. It is the deck's main tool, and it is why
[`qc`](qc.md) buffers a batch of records and then counts:

- **GC and N content** — `core.gcCount` and `core.nCount`, four and two native
  scans per batch
- **adapters** — one scan per probe per batch
- **the quality histogram** — one scan per Phred score, stopping as soon as
  every character in the batch is accounted for, so typical Illumina data costs
  about forty scans rather than the full printable-ASCII range

The results are exact, not sampled or approximate. The test suite checks them
against an independent reference implementation.

## What is still per-base, and why

Three things genuinely need a step per base, and each is opt-in or bounded:

| operation | why it cannot batch |
| --- | --- |
| `qc` with `perBase: true` | position *i* of every read is not a character count |
| `trim.qualityKeepLength` | a sliding window is inherently positional |
| `fastq.qualities` / `meanQuality` / `minQuality` | decodes each score individually |

For per-read work this is fine — 150 steps per read is nothing. It only matters
when it becomes 150 steps × 150 million reads.

## Budgeting a real run

Take **~16 100 reads/second** for a full QC pass as the working figure. A
30× human WGS run of roughly 150 M reads is therefore around three hours of QC
if you read all of it — and there is no reason to.

Subsample. QC statistics are stable long before the file is:

```jennifer
def opts as qc.Options init qc.options();
$opts.everyNth = 20;         # spread across the whole file
$opts.maxReads = 200000;     # and stop at 200k
```

That is a few seconds, and it is what FastQC does too.

For trimming there is no equivalent shortcut — every read has to be written —
so a full trimming pass over a large run is genuinely long. That is the point at
which [`pipeline`](pipelines.md) should be calling `fastp` and this deck should
be parsing its report.

## Memory is bounded by the batch, not the file

Streaming keeps the *file* out of memory, but the QC batch is still real memory,
and it is sized in records. On short reads that is nothing; on long reads it is
not:

| input | batch bound | peak RSS |
| --- | --- | --- |
| 50 000 × 150 bp | 2000 records (300 KB) | 34 MB |
| 3 000 × 30 kb, records only | 2000 records (120 MB) | 657 MB |
| 3 000 × 30 kb, records + bases | 4 MB of bases | 74 MB |

`qc.Options.maxBatchBases` is the second bound, and it is why the third row
exists. Nothing about the results changes — the same file produces the same
statistics either way.

The lesson generalises: when you buffer records to hand a batch to Go, bound the
buffer in **bytes**, not records. Read lengths across NGS span three orders of
magnitude.

## Struct copies in hot loops

Jennifer has value semantics: assignment and argument passing copy. In a module,
a parameter the body never writes is passed by borrow instead, so read-only
arguments are free — but a struct you *build up* per record is not.

An early version of `qc.analyze` accumulated into the `Report` struct itself,
which deep-copied a 94-element histogram and two maps **per record**. Moving the
hot counters to plain locals and assembling the `Report` at the end took the
default path from 6.5 s to 4.0 s for the same work.

The rule that follows: in a per-record loop, accumulate into scalars and
locals. Batch-level state can live in a struct, because it is copied per batch —
2 000× less often.

## Building lists

Use the append sugar:

```jennifer
$xs[] = item;                      # amortised O(n) for n appends
$xs = lists.push($xs, item);       # O(n²) - copies the whole list every pass
```

Two related traps, both of which the deck hit:

- **A map element cannot be appended to in place.** `$m[k][] = v` is a parse
  error, and the read-modify-write workaround is quadratic. `intervals.byChromosome`
  sorts and groups runs instead, in one pass.
- **`binary.concat` in a loop** is O(n²) for the same reason. The deck does not
  accumulate bytes; it accumulates a `list of string` and joins once.

## Parallelism

Two kinds are available, and they solve different problems.

`pipeline.runParallel` fans **external processes** across cores. This is the one
that matters for throughput, because the work is happening in someone else's
optimised C.

`spawn` runs **Jennifer** concurrently. It deep-copies its enclosing scope at
launch, so there are no shared-memory races — and no shared accumulators either.
It suits per-sample work that returns a small summary; it does not turn a
per-base loop into a fast one.

## What the profiler taught this deck

`jennifer profile` reports hits and wall-clock per source position, and
`--allocs` reports every deep copy. Both were run against 200 000-record
fixtures while tuning this deck, and four findings did most of the work. They
generalise to any Jennifer code with a hot loop.

### Binding a list element copies the struct

This is the one that costs the most and looks the most innocent:

```jennifer
for (def i in 1..len($ivs)) {
    def current as Interval init $ivs[$i];   # deep-copies the struct, every pass
    if ($current.start < $limit) { ... }
}
```

The allocation profile showed 400 000 eager copies from a single line like that
inside `intervals.merge`. Reading fields off the indexed element instead has the
same meaning and allocates nothing:

```jennifer
for (def i in 1..len($ivs)) {
    if ($ivs[$i].start < $limit) { ... }
}
```

The same applies to `for (def iv in $xs)` over a list of structs: the loop
variable is a copy per element. It is fine over a list of ints, and fine when
the list is small — it is the per-record path where it matters.

### One exception frame per record, not per field

Guarding each numeric column with its own `try` cost about 2 µs per column.
Measured against a bare conversion: `convert.toInt` 1.4 µs, wrapped in a
function 2.5 µs, wrapped in a function with a `try` 3.5 µs.

SAM has five numeric columns, so per-column guarding cost five frames per
record. Wrapping the whole record parse in one `try` and re-raising through
`core.reraise` costs one — and produces a *better* message, because it names the
record rather than just the column.

### A helper call is not free

`core.field($f, 3, "")` is a clean way to read an optional column, and at two
calls per record it was measurably more expensive than an inline
`if (len($f) > 3)`. Helpers earn their keep everywhere except the innermost
loop; that is where you inline them and say why in a comment.

### Check before you sort

`intervals.sort` builds one zero-padded key string per element and hands them to
`lists.sortBy`. Real BED and GFF files usually arrive coordinate-sorted, and
every set operation sorts its inputs, so most of that work was redundant.
A linear scan comparing fields directly — no keys, no allocation — turns the
common case into one pass, and costs nothing when the data *is* unsorted because
it stops at the first pair out of order.

Measured back to back on a 200 000-feature BED, read-merge-measure takes 16.9 s
on shuffled input and 12.5 s on the same features coordinate-sorted — and the
sorted case allocates 480 MB less, because the sort's output list never exists.

## Where to spend effort

In order:

1. **Do not read what you do not need.** Subsample, or filter early in the
   stream.
2. **Count in Go.** `core.countIn`, `strings.indexOf`, `strings.split` — anything
   that pushes the loop below the interpreter.
3. **Do not copy in the hot loop.** Index list elements rather than binding them;
   keep per-record accumulators in scalars, not structs.
4. **Profile before optimising, and measure after.** Every change above came
   from `jennifer profile`. Three changes that looked obviously good beforehand
   did not survive measurement: a cheaper sort key bought 9%, building the key
   with `strings.join` was *slower* than plain concatenation, and a
   single-character pre-filter in the adapter matcher — which skips a
   slice-and-compare three times in four — came out 18% slower than the
   straightforward loop it replaced, because the slice it added cost as much as
   the one it avoided. All three were reverted. In an interpreter where every
   operation costs about the same, "does less work" and "runs fewer statements"
   are different claims, and only the second one is reliably faster.
5. **Shell out.** If the work is genuinely per-base over a whole run, that is
   what the native tools are for, and [driving them](pipelines.md) is what this
   deck is best at.
