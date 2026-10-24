# Coordinate systems

Read this before writing interval code. Mixing the two conventions shifts every
downstream result by one base and never raises an error.

## The two conventions

| format | convention | `chr1:100-200` means |
| --- | --- | --- |
| BED, `intervals.Interval` | 0-based, half-open `[start, end)` | bases 101-200 in 1-based terms, 100 bases |
| GFF3, GTF | 1-based, inclusive | bases 100-200, 101 bases |
| SAM `POS` | 1-based | the leftmost aligned base |
| VCF `POS` | 1-based | the first base of `REF` |

A BED interval's length is `end - start`; a GFF feature's is `end - start + 1`.

## Converting

Records keep their file's coordinates: `gff.Feature.start` is the number in
column 4, `sam.Record.pos` the number in column 4 of the SAM line. Nothing is
normalised on read, so a record written back matches what was read.

Conversion is an explicit call, and these calls are the way into the interval
engine:

```jennifer
def a as intervals.Interval init gff.toInterval($feature);   # start - 1
def b as intervals.Interval init sam.toInterval($alignment); # pos - 1, end from CIGAR
def c as intervals.Interval init vcf.toInterval($variant);   # pos - 1, spans REF
def f as gff.Feature init gff.fromInterval($iv, "ngs", "region");  # start + 1
```

`gff.toInterval` names the interval with `gff.name`; `gff.fromInterval` sets
the `ID` attribute from the interval's name and leaves `phase` at -1.
`vcf.toInterval` names it with the variant ID, and a missing QUAL becomes a
score of `0.0`. BED needs no conversion: `bed.parse` returns an
`intervals.Interval`.

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

Only the start shifts: a 1-based inclusive end and a 0-based exclusive end are
the same number. Subtracting one from both is a common wrong fix.

## What `sam.toInterval` does with CIGAR

Deletions (`D`) and skipped regions (`N`) consume reference; insertions, soft
and hard clips and padding do not.

```jennifer
# 10M5D10M starting at POS 101
def iv as intervals.Interval init sam.toInterval($rec);
$iv.start;   # 100
$iv.end;     # 125  - 25 reference bases, not the 20 query bases
```

When the CIGAR is `*`, the span falls back to the sequence length. The interval
is named after the read, scored with MAPQ and stranded from FLAG `0x10`. An
unmapped record (FLAG `0x4`, or RNAME `*`) raises `ngs.sam`; filter with
`sam.isUnmapped` first.

## Sort order

`intervals.sort` orders by chromosome as text, then by start, then by end, the
order of `sort -k1,1 -k2,2n` and bedtools: `chr10` sorts before `chr2`. A set
sorted here can go to those tools and back without a re-sort. For karyotypic
order, sort afterwards.

The order comes from `intervals.sortKey`, which zero-pads coordinates to
`intervals.KEY_WIDTH` (12) digits. A coordinate of `intervals.MAX_COORDINATE`
(10^12) or more would not fit, so `intervals.make` and `intervals.named` reject
it with `ngs.intervals`.

## Book-ended intervals do not overlap

`[10, 20)` and `[20, 30)` touch but share no base: `intervals.overlaps` returns
`false` and `intervals.overlapLength` returns `0`.

`intervals.merge` with a `maxGap` of 0 joins them anyway, as
`bedtools merge -d 0` does. A larger `maxGap` bridges real gaps.
