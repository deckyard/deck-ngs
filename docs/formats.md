# Reading and writing formats

Every reader is opened the same way, streams, decompresses gzip transparently
and can read standard input.

| module | format | record type | coordinates as stored |
| --- | --- | --- | --- |
| `fastq` | FASTQ | `fastq.Read` | n/a |
| `bed` | BED | `intervals.Interval` | 0-based, half-open |
| `gff` | GFF3 and GTF | `gff.Feature` | 1-based, inclusive |
| `sam` | SAM | `sam.Record` | 1-based (`pos`) |
| `vcf` | VCF | `vcf.Variant` | 1-based (`pos`) |

Records keep their file's coordinates; [Coordinate systems](coordinates.md)
covers converting them.

## Two cursor shapes

FASTQ has no comment lines and every record is four lines, so its reader uses
`hasNext` / `next`:

```jennifer
def r as fastq.Reader init fastq.open("reads.fastq.gz");
while (fastq.hasNext($r)) {
    def rec as fastq.Read init fastq.next($r);
    io.printf("%s %d\n", $rec.id, fastq.length($rec));
}
fastq.closeStrict($r);
```

The tab-delimited formats can carry comments, headers and directives between
records: a `track` line in BED, a `###` directive or a trailing `##FASTA`
section in GFF3. Skipping them can consume the last line of the file, and a
reader cannot un-read a line ([why](streaming.md#value-semantics-and-handles)),
so a `hasNext` peek cannot answer correctly. These readers return a `Next` enum
instead:

```jennifer
def r as bed.Reader init bed.open("peaks.bed");
def going as bool init true;
while ($going) {
    def step as bed.Next init bed.next($r);
    match ($step) {
        when Feature(f) { io.printf("%s\n", intervals.toBed3($f.value)); }
        when Done { $going = false; }
    }
}
bed.closeStrict($r);
```

The `match` subject must be a variable of the enum type; `match (bed.next($r))`
is rejected. The record variant is `bed.Next.Feature`, `gff.Next.Record`,
`sam.Next.Alignment` or `vcf.Next.Record`; the record is in its `value` field.

## FASTQ

`fastq.Read` holds `id`, `desc`, `seq` and `qual`. The header is split at the
first space: `@r1 length=150` gives `id` `r1` and `desc` `length=150`.

```jennifer
def rec as fastq.Read init fastq.parse("@r1 note", "ACGTN", "+", "IIII#");

fastq.length($rec);                       # 5
fastq.gcPercent($rec);                    # 40.0
fastq.ambiguousCount($rec);               # 1
fastq.meanQuality($rec, fastq.PHRED33);   # arithmetic mean Phred, as FastQC and fastp
fastq.minQuality($rec, fastq.PHRED33);    # -1 for an empty read
fastq.qualities($rec, fastq.PHRED33);     # list of int, one per base

fastq.subseq($rec, 1, 4);                 # 0-based, half-open; slices seq and qual
fastq.reverseComplement($rec);            # reverses qual to match
fastq.format($rec);                       # four lines, newline-terminated
fastq.toFasta($rec);                      # two-line FASTA
```

`fastq.PHRED33` (33) is the Sanger / Illumina 1.8+ offset; `fastq.PHRED64` (64)
is legacy Illumina 1.3-1.7. `qualities`, `meanQuality` and `minQuality` decode
one base per interpreter step; for whole-file statistics use [`qc`](qc.md).

`parse` and `next` raise `ngs.fastq` on a header not starting with `@`, a
separator not starting with `+`, a quality byte outside `!`..`~` (ASCII
33..126), a sequence and quality of different byte lengths, or a truncated
record; `next` past the end raises `ngs.fastq` too. Blank lines are skipped:
`hasNext` reads ahead to the next non-blank line and keeps it for `next`, so a
file ending in empty lines ends cleanly.

`format` always writes the separator as a bare `+`. It raises `ngs.fastq` on a
record that would not read back as itself: whitespace in `id`, a line break in
`desc` (a tab is kept, as `samtools fastq -T` writes tags there), a tab or line
break in `seq`, a bad quality byte or a length mismatch.

[Recipes](recipes.md#convert-fastq-to-fasta) has a streaming FASTQ-to-FASTA
conversion.

## BED

BED records are `intervals.Interval` values, so they go straight into the
[interval engine](intervals.md). The first six columns are read; thickStart,
blocks and later columns are ignored. A missing name reads as `""` and a missing
strand as `"."`; `bed.format` writes both back as `.`, and a whole-number score
without a decimal part.

```jennifer
def iv as intervals.Interval init bed.parse("chr1\t100\t200\tpeak1\t55\t-");
bed.format($iv);                # chr1  100  200  peak1  55  -
```

Annotation is usually small next to sequencing data, so the whole-file calls
are the common path:

```jennifer
def peaks as list of intervals.Interval init bed.readAll("peaks.bed");
bed.writeAll("out.bed.gz", $peaks);    # gzipped because of the name
```

Readers skip blank lines, `#` comments and `track` / `browser` lines. A record
with fewer than three columns, or a non-numeric coordinate, raises `ngs.bed`;
negative or inverted coordinates, or a strand other than `+`, `-` or `.`, raise
`ngs.intervals`. `bed.write` appends one interval to an open sink. `bed.format`
raises `ngs.bed` when a chromosome, name or strand contains a tab or line
break, and keeps a score of 10^15 or more in float notation.

## GFF3 and GTF

One module reads both. They share nine tab-delimited columns and differ in the
ninth: GFF3 writes `key=value;key=value` with percent-escaping, GTF writes
`key "value"; key "value";`. `gff.parse` detects the style per record from the
first pair. A GTF value written unquoted, as GENCODE does for `level 2`, is
still read as GTF.

```jennifer
def f as gff.Feature init gff.parse($line);

$f.seqid; $f.origin; $f.kind;      # the seqid, source and type columns
$f.start; $f.end;                  # 1-based, inclusive, as written
$f.strand; $f.phase;               # phase is -1 when absent
$f.score;                          # 0.0 when absent

gff.attribute($f, "gene_id", "");
gff.identifier($f);                # ID, else transcript_id, else gene_id
gff.name($f);                      # Name, else gene_name, else identifier
```

GFF3 values are percent-decoded on read. On write, only the characters GFF3
reserves are escaped: tab, newline, carriage return and `%`, plus `;`, `=` and
`&` in column 9. `,` is not escaped: it separates the values of a multi-valued
attribute, so `Parent=t1,t2` round-trips as two parents. Spaces stay as they
are. An empty attribute map, or a ninth column of `.`, is written as `.`. A
score of `0.0` is written as `.`.

A multi-valued attribute is one string joined with `,`. A key that appears
more than once - GTF's `tag "basic"; tag "CCDS";`, a GFF3 `Dbxref` written
twice - reads as its values joined with `,` in file order (`"basic,CCDS"`); a
GTF write emits that as one quoted value, not as repeated keys. A decoded
`%2C` cannot be told from the separator. GTF attributes are split on `;`
outside double quotes, so `note "a; b"` stays one value.

`gff.format` raises `ngs.gff` on a tab or line break in columns 1-8, an
attribute key that is empty or contains `;`, `=`, `"` or whitespace, and a GTF
value containing `"`, a tab or a line break, which GTF has no escape for.
`gff.toInterval` maps the GFF3 strand `?` (relevant but unknown) to `.`.

```jennifer
def exons as list of intervals.Interval init gff.readIntervals("genes.gff3", "exon");
def all as list of intervals.Interval init gff.readIntervals("genes.gff3", "");
```

Readers stop at a `##FASTA` section. A record with fewer than nine columns
raises `ngs.gff`. `gff.readAll` closes its reader with `close`, not
`closeStrict`, so a truncated `.gz` annotation is not reported.

## SAM

SAM only; read BAM or CRAM [through `samtools view -h`](recipes.md#read-bam-through-samtools).

```jennifer
def rec as sam.Record init sam.parse($line);

sam.isUnmapped($rec);
sam.isReverse($rec);
sam.isSecondary($rec);
sam.isSupplementary($rec);
sam.isDuplicate($rec);
sam.isPrimary($rec);          # mapped, not secondary, not supplementary
sam.strandOf($rec);           # "+", "-", or "." when unmapped
sam.hasFlag($rec, sam.FLAG_PROPER_PAIR);   # true when any bit of the mask is set

sam.tagValue($rec, "NM", "");  # value with its TAG:TYPE: prefix stripped
```

The flag constants, `sam.FLAG_PAIRED` (`0x1`) to `sam.FLAG_SUPPLEMENTARY`
(`0x800`), are the SAM specification's FLAG bits. A record needs at
least eleven columns; columns after the eleventh are kept verbatim in `tags`.

`sam.parseCigar` splits a CIGAR into `sam.CigarOp` values and raises `ngs.sam`
on a malformed string. The two span functions:

```jennifer
sam.referenceLength("10M5D10M");   # 25 - counts M D N = X
sam.queryLength("10M2I5S");        # 17 - counts M I S = X
```

`sam.toInterval` derives the end from `referenceLength`, so a spliced or
deleted alignment gets its reference span rather than its read length
([details](coordinates.md#what-samtointerval-does-with-cigar)). A CIGAR of `*`,
or one that consumes no reference bases, gives a one-base interval at POS, as
htslib's `bam_endpos` does. `sam.format` raises `ngs.sam` when a field or tag
contains a tab or line break.

### The header

`sam.open` reads the `@` header block with a separate short pass over the
front of the file and keeps it on `Reader.header`; the record stream then
starts at byte zero and skips the `@` lines. `sam.referenceSizes` turns the
`@SQ` lines into the chromosome-size map
[`intervals.complement`](intervals.md#complement) takes:

```jennifer
def r as sam.Reader init sam.open("aln.sam");
def sizes as map of string to int init sam.referenceSizes($r);
```

On standard input `header` is empty: a stream cannot be re-read, and reading
the header first would consume the first alignment. `sam.headerOf(path)` reads
only the header of a file; it closes with `close`, since stopping early makes a
gzip decompressor exit on `SIGPIPE`.

## VCF

```jennifer
def v as vcf.Variant init vcf.parse($line);

$v.chrom; $v.pos; $v.id; $v.ref;
$v.alt;                        # list of string, empty when "."
$v.qual;                       # -1.0 when the column is "."
$v.filter;                     # verbatim

vcf.passed($v);                # FILTER is PASS, "." or empty
vcf.isSnv($v);                 # REF and every ALT one of A C G T N (upper case)
vcf.isIndel($v);               # an ALT of different length, symbolic and breakend ALTs excluded
vcf.info($v, "DP", "0");
vcf.hasInfo($v, "DB");         # a bare flag maps to ""
```

`vcf.parseInfo` and `vcf.formatInfo` convert the INFO column to and from a
map; a value runs from the first `=` to the end of the entry. Values are
percent-encoded on write and decoded on read as VCF 4.3 specifies (`%` `:` `;`
`=` CR LF tab; `%2C` is decoded too), so a value may contain `;`. `,` is not
encoded: a multi-valued field such as `AF=0.1,0.2` is one string. A key that
appears twice raises `ngs.vcf`, as does writing an empty key or one containing
`;`, `=` or whitespace.

Per-sample values are read by index into the reader's `samples` list:

```jennifer
def r as vcf.Reader init vcf.open("calls.vcf.gz");
for (def i in 0..len($r.samples)) {
    io.printf("%s %s\n", $r.samples[$i], vcf.genotype($v, $i));   # "./." when absent
}
vcf.sampleValue($v, 0, "DP", ".");
```

`vcf.open` captures `meta` (the `##` lines) and `samples` (from the `#CHROM`
line) with a separate pass, as `sam.open` does; on standard input both are
empty. A record with fewer than eight columns raises `ngs.vcf`.

`vcf.format` reproduces a parsed record, minimal eight-column records included,
with one exception: QUAL is re-rendered from a float, so an integer QUAL such
as `60` is written as `60.0`. An empty FILTER is written `.`, and an empty
FORMAT is written `.` when the record has sample columns, so they stay in
place. A tab or line break in any other column raises `ngs.vcf`.
