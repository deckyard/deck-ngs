# Coordinate systems

This is the classic source of silent, off-by-one bioinformatics bugs. The deck
is explicit about it, and this page is the one to read before writing interval
code.

## The two conventions

| format | convention | `chr1:100-200` means |
| --- | --- | --- |
| BED, `intervals.Interval` | **0-based, half-open** `[start, end)` | bases 101–200 in 1-based terms, 100 bases |
| GFF3, GTF | **1-based, inclusive** | bases 100–200, 101 bases |
| SAM `POS` | **1-based** | the leftmost aligned base |
| VCF `POS` | **1-based** | the first base of `REF` |

A BED interval's length is `end - start`. A GFF feature's length is
`end - start + 1`. Getting these confused shifts every downstream result by one
base and never raises an error.

## The rule

**Records keep their file's native coordinates.** `gff.Feature.start` is exactly
the number in column 4; `sam.Record.pos` is exactly the number in column 4 of
the SAM record. Nothing is silently normalised on read, so a record you print
back out matches what you read in.

**Converting is always an explicit call**, and those calls are the only correct
way to move a record into the interval engine:

```jennifer
def a as intervals.Interval init gff.toInterval($feature);   # start - 1
def b as intervals.Interval init sam.toInterval($alignment); # pos - 1, end from CIGAR
def c as intervals.Interval init vcf.toInterval($variant);   # pos - 1, spans REF
```

Going back the other way:

```jennifer
def f as gff.Feature init gff.fromInterval($iv, "ngs", "region");  # start + 1
```

BED needs no conversion in either direction — `bed.parse` returns an
`intervals.Interval` directly, because the conventions already match.

## Worked example

A GFF exon written as `101 200`:

```jennifer
def f as gff.Feature init gff.parse("chr1\tsrc\texon\t101\t200\t.\t+\t.\tID=e1");
$f.start;                        # 101  - as written in the file
$f.end;                          # 200

def iv as intervals.Interval init gff.toInterval($f);
$iv.start;                       # 100  - converted
$iv.end;                         # 200  - unchanged
intervals.length($iv);           # 100  - both conventions agree on the count
```

The end coordinate does not move. A 1-based inclusive end and a 0-based
exclusive end are the same number — only the start shifts. That is why
converting by subtracting one from *both* is a common and wrong fix.

## What `sam.toInterval` does with CIGAR

An alignment's reference span is not its read length: deletions and skipped
regions consume reference, insertions and clips do not.

```jennifer
# 10M5D10M starting at POS 101
def iv as intervals.Interval init sam.toInterval($rec);
$iv.start;   # 100
$iv.end;     # 125  - 25 reference bases, not the 20 query bases
```

`sam.referenceLength` counts `M D N = X`; `sam.queryLength` counts `M I S = X`.
When the CIGAR is `*`, `toInterval` falls back to the sequence length.

An unmapped record has no place on the reference, so `sam.toInterval` raises
`ngs.sam` rather than inventing one. Filter with `sam.isUnmapped` first.

## Sort order

`intervals.sort` orders **lexicographically by chromosome, then numerically by
start, then by end**. That matches `sort -k1,1 -k2,2n` and bedtools' default
order, which means `chr10` sorts before `chr2`.

This is deliberate: it is what the surrounding tools do, and an interval set
sorted here can be handed to them and back without a resort. If you need
karyotypic order, sort your own way after the fact.

Internally the ordering comes from `intervals.sortKey`, which builds a
lexicographically sortable string with the coordinates zero-padded to
`intervals.KEY_WIDTH` digits. It is exported because `lists.sortBy` needs it as a
`func` value; you rarely call it yourself.

## Book-ended intervals do not overlap

With half-open coordinates, `[10, 20)` and `[20, 30)` touch but share no base.
`intervals.overlaps` returns `false` for them and `intervals.overlapLength`
returns `0`.

`intervals.merge` is the exception, and matches `bedtools merge -d 0`: it joins
overlapping *and* book-ended features, because a gap of zero is not a gap. Pass
a larger `maxGap` to bridge real gaps.
