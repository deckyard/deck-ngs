# Interval arithmetic

`intervals` implements the common bedtools operations over
`intervals.Interval`, which is 0-based and half-open
([Coordinate systems](coordinates.md) shows how to convert other formats).

```jennifer
def struct Interval {
    chrom as string,
    start as int,      # 0-based, inclusive
    end as int,        # 0-based, exclusive
    name as string,    # "" when unnamed
    score as float,    # 0.0 when absent
    strand as string   # "+", "-", or "."
};
```

## Building and inspecting

```jennifer
def a as intervals.Interval init intervals.make("chr1", 100, 200);
def b as intervals.Interval init intervals.named("chr1", 150, 250, "peak", 55.0, "+");

intervals.length($a);            # 100
intervals.overlaps($a, $b);      # true
intervals.overlapLength($a, $b); # 50
intervals.contains($a, $b);      # false
intervals.intersection($a, $b);  # chr1:150-200, with $a's name, score and strand
intervals.toBed3($a);            # "chr1\t100\t200"
```

`make` builds an unnamed interval with strand `"."`. `make` and `named` raise
`ngs.intervals` on a negative coordinate, an inverted interval (`end < start`)
or an end of `intervals.MAX_COORDINATE` (10^12) or more.
`intervals.intersection` raises `ngs.intervals` when the two do not overlap.

## Set operations

Each takes a query set `a` and a target set `b`. Results follow
`intervals.sort(a)`, not the input order.

| call | bedtools equivalent | returns |
| --- | --- | --- |
| `intersect(a, b)` | `intersect -a A -b B` | one interval per overlapping pair, covering the shared bases, carrying the query's name |
| `countOverlaps(a, b)` | `intersect -c` | `list of int`, one count per query |
| `coverage(a, b)` | covered-bases column of `coverage` | `list of int`, bases of each query covered by `b` |
| `subtract(a, b)` | `subtract` | the pieces of each query not covered by `b` |

```jennifer
def peaks as list of intervals.Interval init bed.readAll("peaks.bed");
def exons as list of intervals.Interval init gff.readIntervals("genes.gff3", "exon");

def shared as list of intervals.Interval init intervals.intersect($peaks, $exons);
def hits as list of int init intervals.countOverlaps($peaks, $exons);
def covered as list of int init intervals.coverage($peaks, $exons);
def outside as list of intervals.Interval init intervals.subtract($peaks, $exons);
```

`coverage` and `subtract` merge the targets first, so a base under two targets
counts once. `intersect` and `countOverlaps` do not: every overlapping pair is
reported.

`subtract` splits a query with a target through its middle into two pieces and
drops a fully covered query. Surviving pieces keep the query's name, score and
strand.

All four sort both sets, then sweep a cursor forward through the targets as the
queries advance; a target ending before the current query starts cannot overlap
a later query. Cost is linear in the input plus the number of reported
overlaps, after the sort, which runs in Go through `lists.sortBy`.
`intervals.sort` returns input that is already in
order without sorting it.

## Merging and measuring

```jennifer
intervals.merge($xs, 0);      # overlapping and book-ended, like bedtools merge -d 0
intervals.merge($xs, 100);    # also joins features up to 100 bp apart
intervals.totalLength($xs);   # bases covered by the union
```

`merge` works per chromosome and returns sorted, disjoint intervals with name
`""`, score `0.0` and strand `"."`: a merged run is not any one of its inputs.

## Complement

The uncovered stretches of each chromosome, like `bedtools complement -g`. The
chromosome sizes can come from a SAM header:

```jennifer
def sizes as map of string to int init sam.referenceSizes($reader);
def gaps as list of intervals.Interval init intervals.complement($covered, $sizes);
```

Only chromosomes named in `sizes` appear in the result, which is sorted.

## Growing and grouping

```jennifer
intervals.slop($iv, 1000, 1000, 0);        # 1 kb each side, clamped at 0 only
intervals.slop($iv, 1000, 1000, $limit);   # also clamped at the chromosome end

intervals.byChromosome($xs);               # map of string to list of Interval
```

`slop` is `bedtools slop`; a `limit` of 0 or less means no upper clamp.
`byChromosome` sorts first, so each group is sorted by start and the grouping
is one pass. For why these functions avoid per-element copies, see
[Performance](performance.md#struct-copies-in-hot-loops).
