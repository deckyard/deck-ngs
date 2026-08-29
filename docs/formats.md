# Reading and writing formats

Five formats, one shape: open a cursor, pull records, close it. Every reader
streams, decompresses gzip transparently, and can read standard input.

| module | format | record type | coordinates |
| --- | --- | --- | --- |
| `fastq` | FASTQ | `fastq.Read` | n/a |
| `bed` | BED | `intervals.Interval` | 0-based, half-open |
| `gff` | GFF3 and GTF | `gff.Feature` | 1-based, inclusive |
| `sam` | SAM | `sam.Record` | 1-based (`pos`) |
| `vcf` | VCF | `vcf.Variant` | 1-based (`pos`) |

Records keep their file's native coordinates. See
[Coordinate systems](coordinates.md) before mixing them.

## Two cursor shapes

FASTQ has no comment lines: every record is exactly four lines, so
`hasNext`/`next` works and reads naturally.

```jennifer
def r as fastq.Reader init fastq.open("reads.fastq.gz");
while (fastq.hasNext($r)) {
    def rec as fastq.Read init fastq.next($r);
    io.printf("%s %d\n", $rec.id, fastq.length($rec));
}
fastq.closeStrict($r);
```

The tab-delimited formats can carry comments, headers and directives *between*
records - a `track` line in BED, a `###` directive or a trailing `##FASTA`
section in GFF3. Skipping those can consume the last line of the file, which a
`hasNext` peek cannot express: it would have to answer "is there another record"
without reading, and a reader cannot un-read a line
([why](streaming.md#value-semantics-and-handles)).

So those cursors report end-of-stream through an enum, which `match` checks
exhaustively:

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

The `match` subject must be a variable of the enum type - `match (bed.next($r))`
is rejected, because the subject has to be statically an enum.

The variant name differs per module, since the payload does: `bed.Next.Feature`,
`gff.Next.Record`, `sam.Next.Alignment`, `vcf.Next.Record`.

## FASTQ

`fastq.Read` is `{id, desc, seq, qual}`. The header is split at the first space:
`@r1 length=150` gives `id` of `r1` and `desc` of `length=150`.

```jennifer
def rec as fastq.Read init fastq.parse("@r1 note", "ACGTN", "+", "IIII#");

fastq.length($rec);                       # 5
fastq.gcPercent($rec);                    # 40.0
fastq.ambiguousCount($rec);               # 1
fastq.meanQuality($rec, fastq.PHRED33);   # arithmetic mean Phred
fastq.minQuality($rec, fastq.PHRED33);
fastq.qualities($rec, fastq.PHRED33);     # list of int, one per base

fastq.subseq($rec, 1, 4);                 # slices seq and qual together
fastq.reverseComplement($rec);            # reverses qual to match
fastq.format($rec);                       # back to four lines
fastq.toFasta($rec);                      # two-line FASTA
```

`qualities` and `meanQuality` decode the quality string, which costs one
interpreter step per base. For whole-file statistics use [`qc`](qc.md), which
counts in Go-speed batches instead.

Parsing validates the shape: a header not starting with `@`, a separator not
starting with `+`, or a sequence and quality of different lengths all raise
`ngs.fastq`.

### Converting as a stream

FASTQ to FASTA, without either file in memory:

```jennifer
def r as fastq.Reader init fastq.open("in.fastq.gz");
def k as source.Sink init source.create("out.fasta.gz");
while (fastq.hasNext($r)) {
    source.writeText($k, fastq.toFasta(fastq.next($r)));
}
fastq.closeStrict($r);
source.closeSink($k);
```

## BED

BED maps directly onto `intervals.Interval` - both are 0-based and half-open -
so a parsed feature drops straight into [the interval engine](intervals.md).
The first six columns are read; anything further (thickStart, blocks) is ignored.

```jennifer
def iv as intervals.Interval init bed.parse("chr1\t100\t200\tpeak1\t55\t-");
bed.format($iv);                # chr1  100  200  peak1  55  -
```

For feature sets that fit in memory - the usual case, since annotation is small
next to sequencing data - skip the cursor entirely:

```jennifer
def peaks as list of intervals.Interval init bed.readAll("peaks.bed");
bed.writeAll("out.bed.gz", $peaks);
```

`readAll` skips blank lines, `#` comments, and `track` / `browser` lines.

## GFF3 and GTF

One module handles both. They share nine tab-delimited columns and differ only
in the ninth: GFF3 writes `key=value;key=value` with percent-escaping, GTF writes
`key "value"; key "value";`. `gff.parse` detects the style per record, so a mixed
directory needs no flag.

```jennifer
def f as gff.Feature init gff.parse($line);

$f.seqid; $f.origin; $f.kind;      # 'source' and 'type' columns
$f.start; $f.end;                  # 1-based, inclusive, as written
$f.strand; $f.phase;               # phase is -1 when absent

gff.attribute($f, "gene_id", "");
gff.identifier($f);                # ID, else transcript_id, else gene_id
gff.name($f);                      # Name, else gene_name, else identifier
```

GFF3 attribute values are percent-decoded on read and re-escaped on write, so
a `Note` containing `;` or `=` survives a round trip.

Reading a whole annotation, optionally filtered by feature type, and converting
to intervals in one call:

```jennifer
def exons as list of intervals.Interval init gff.readIntervals("genes.gff3", "exon");
def all as list of intervals.Interval init gff.readIntervals("genes.gff3", "");
```

The cursor stops at a `##FASTA` section rather than trying to parse the sequence
block that follows it.

## SAM

SAM only. BAM and CRAM are [out of scope](index.md#out-of-scope) - pipe
`samtools view -h` in instead:

```sh
samtools view -h aln.bam | jennifer run myscript.j -
```

```jennifer
def rec as sam.Record init sam.parse($line);

sam.isUnmapped($rec);
sam.isReverse($rec);
sam.isSecondary($rec);
sam.isSupplementary($rec);
sam.isDuplicate($rec);
sam.isPrimary($rec);          # mapped, not secondary, not supplementary
sam.strandOf($rec);           # "+", "-", or "." when unmapped
sam.hasFlag($rec, sam.FLAG_PROPER_PAIR);

sam.tagValue($rec, "NM", "");  # optional tag, TAG:TYPE: prefix stripped
```

CIGAR strings are parsed into `sam.CigarOp{length, op}`, with the two span
functions that matter:

```jennifer
sam.referenceLength("10M5D10M");   # 25 - counts M D N = X
sam.queryLength("10M2I5S");        # 17 - counts M I S = X
```

`sam.toInterval` uses `referenceLength` to derive the end position, so a spliced
or deleted alignment gets its true reference span rather than the read length.

### The header

`sam.open` captures the `@` header block on the `Reader`, and
`sam.referenceSizes` turns the `@SQ` lines into the chromosome-size map that
[`intervals.complement`](intervals.md#complement) wants:

```jennifer
def r as sam.Reader init sam.open("aln.sam");
def sizes as map of string to int init sam.referenceSizes($r);
```

The header is read by a short separate pass over the front of the file, so the
record stream still starts clean. On **stdin** `header` is empty - a stream
cannot be re-read, and reading it eagerly would consume the first alignment with
no way to put it back. `sam.headerOf(path)` reads a header on its own.

## VCF

```jennifer
def v as vcf.Variant init vcf.parse($line);

$v.chrom; $v.pos; $v.id; $v.ref; $v.alt;   # alt is a list
$v.qual;                                    # -1.0 when the column is "."
$v.filter;

vcf.passed($v);                # PASS or "." both count as passing
vcf.isSnv($v);                 # ref and every alt a single base
vcf.isIndel($v);               # a length difference, symbolic alleles excluded
vcf.info($v, "DP", "0");
vcf.hasInfo($v, "DB");         # a bare flag maps to ""
```

Per-sample values are read by index into the reader's `samples` list:

```jennifer
def r as vcf.Reader init vcf.open("calls.vcf.gz");
for (def i in 0..len($r.samples)) {
    io.printf("%s %s\n", $r.samples[$i], vcf.genotype($v, $i));
}
vcf.sampleValue($v, 0, "DP", ".");
```

`vcf.open` captures `meta` (the `##` lines) and `samples` (from `#CHROM`) the
same way `sam.open` does, with the same stdin caveat.

`vcf.format` round-trips a parsed record byte-for-byte for the common shapes,
including minimal eight-column records with no genotypes.
