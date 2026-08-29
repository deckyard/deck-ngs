# Interval arithmetic

The `intervals` module is the bedtools-lite layer: overlap, merge, intersect,
subtract, coverage and complement over genomic features.

All of it works on `intervals.Interval`, which is **0-based and half-open** —
see [Coordinate systems](coordinates.md) for how to get records from other
formats into it.

```jennifer
def struct Interval {
    chrom as string,
    start as int,      # 0-based, inclusive
    end as int,        # 0-based, exclusive
    name as string,    # "" when unnamed
    score as float,
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
intervals.intersection($a, $b);  # chr1:150-200, named after $a
intervals.toBed3($a);            # "chr1\t100\t200"
```

`make` and `named` validate: a negative coordinate or an inverted interval
(`end < start`) raises `ngs.intervals` rather than producing something that
silently misbehaves later.

## Set operations

Each of these takes a query set `a` and a target set `b`, and returns results
aligned with `intervals.sort(a)` — sorted, not in input order.

| call | bedtools equivalent | returns |
| --- | --- | --- |
| `intersect(a, b)` | `intersect -a A -b B` | one interval per overlapping pair, covering the shared bases |
| `countOverlaps(a, b)` | `intersect -c` | `list of int`, one count per query |
| `coverage(a, b)` | the covered-bases column of `coverage` | `list of int`, bases of each query covered by `b` |
| `subtract(a, b)` | `subtract` | the pieces of each query not covered by `b` |

```jennifer
def peaks as list of intervals.Interval init bed.readAll("peaks.bed");
def exons as list of intervals.Interval init gff.readIntervals("genes.gff3", "exon");

def shared as list of intervals.Interval init intervals.intersect($peaks, $exons);
def hits as list of int init intervals.countOverlaps($peaks, $exons);
def covered as list of int init intervals.coverage($peaks, $exons);
def outside as list of intervals.Interval init intervals.subtract($peaks, $exons);
```

`coverage` merges the target set first, so a base covered by two overlapping
targets counts once. `intersect` does not — it reports every overlapping pair,
which is what you want when the pairing itself is the information.

`subtract` splits: a query with a target through its middle yields two pieces,
and a fully covered query yields none. Names, scores and strands are carried
onto the surviving pieces.

### Why these are fast

All four run a **sorted sweep**, not a nested scan. Both sets are sorted once,
then a cursor walks the targets forward as the queries advance — safe because a
target ending before the current query's start cannot overlap any later query
either. That makes them linear in the input plus the number of reported
overlaps, rather than quadratic.

Sorting itself happens in Go via `lists.sortBy`, so it is not the bottleneck.

## Merging and measuring

```jennifer
intervals.merge($xs, 0);      # overlapping and book-ended, like bedtools merge -d 0
intervals.merge($xs, 100);    # also joins features up to 100 bp apart
intervals.totalLength($xs);   # bases covered by the union, counting overlaps once
```

`merge` returns sorted, disjoint intervals with names cleared — a merged feature
is not any one of its inputs, so keeping one input's name would be a lie.

## Complement

The gaps between features, across whole chromosomes. Needs a chromosome-size
map, which you can get from a SAM header:

```jennifer
def sizes as map of string to int init sam.referenceSizes($reader);
def gaps as list of intervals.Interval init intervals.complement($covered, $sizes);
```

Only chromosomes named in `sizes` appear in the result.

## Growing and grouping

```jennifer
intervals.slop($iv, 1000, 1000, 0);        # 1 kb each side, unclamped
intervals.slop($iv, 1000, 1000, $limit);   # clamped at 0 and the chromosome end

intervals.byChromosome($xs);               # map of string to list of Interval
```

`byChromosome` sorts first and groups the runs, so each group is itself sorted
and the whole thing costs one pass. Building the same map by appending per
element would be quadratic, because a map element cannot be appended to in
place.

## A note on value semantics

These functions take lists and return new ones; nothing mutates its argument. In
a module, a parameter the body never writes is passed by borrow rather than
deep-copied, so handing a large interval list to `intersect` does not copy it.

Building lists is the case to watch. Use the `$xs[] = item;` append sugar, which
is amortised linear; `$xs = lists.push($xs, item)` in a loop copies the whole
list every pass and turns an O(n) build into O(n²).
