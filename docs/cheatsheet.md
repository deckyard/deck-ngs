# Cheatsheet

Every export, by module. Jennifer does not declare return types, so the
`returns` column gives them. Every error is an `Error` whose `kind` is
`ngs.<module>` unless noted.

## `core`

Errors, field parsing, Go-speed counting and shell quoting. Every other module
imports it.

| constant | value | notes |
| --- | --- | --- |
| `core.ERROR_KIND` | `"ngs"` | prefix of every error kind the deck raises |
| `core.TAB` | `"\t"` | field separator of every text format here |

| function | returns | notes |
| --- | --- | --- |
| `core.fail(what, message)` | - | always raises kind `ngs.<what>` |
| `core.countIn(haystack, needle)` | `int` | non-overlapping occurrences, one native scan; 0 for an empty needle |
| `core.parseInt(what, label, raw)` | `int` | raises `ngs.<what>` naming `label` when `raw` is not an integer |
| `core.parseFloat(what, label, raw)` | `float` | raises `ngs.<what>` when `raw` is not a number |
| `core.reraise(what, line, e)` | - | re-raises `e` untouched if its kind starts `ngs.`, else raises `ngs.<what>` naming `line` |
| `core.splitTab(line)` | `list of string` | split on tabs |
| `core.field(fields, i, fallback)` | `string` | column `i` (0-based), or `fallback` when out of range |
| `core.intField(what, fields, i, fallback)` | `int` | `fallback` for a missing, empty or `.` column; raises `ngs.<what>` on a non-integer |
| `core.floatField(what, fields, i, fallback)` | `float` | as `intField`, for floats |
| `core.isSkippable(line, marker)` | `bool` | true for a blank line or one starting with `marker` |
| `core.padNum(n, width)` | `string` | zero-padded to `width`; negative values unpadded |
| `core.shellQuote(s)` | `string` | single-quoted for POSIX `sh`, embedded `'` escaped |
| `core.reverseComplement(seq)` | `string` | ACGTUN and IUPAC codes, case kept; other characters unchanged |
| `core.complementBase(b)` | `string` | complement of one base, or `b` when not a nucleotide code |
| `core.gcCount(seq)` | `int` | G, C, g and c, counted natively |
| `core.nCount(seq)` | `int` | N and n, counted natively |
| `core.hasTool(name)` | `bool` | a regular file of that name on `PATH` (executable bit not checked); walks `PATH` with `fs`, no `exec` needed; a name with a separator is checked as a path |
| `core.sameFile(a, b)` | `bool` | both exist and resolve, symlinks followed, to the same file |

## `source`

Streaming line sources and sinks, gzip in both directions
([Streaming large files](streaming.md)). Errors are `ngs.source`.

| constant | value | notes |
| --- | --- | --- |
| `source.GZIP_MAGIC1` | `0x1f` | first gzip magic byte (RFC 1952) |
| `source.GZIP_MAGIC2` | `0x8b` | second gzip magic byte |

### `source.Source`

An open line source. Copies share the stream: `file` and `proc` are registry
handles.

| field | notes |
| --- | --- |
| `file as fs.File` | the handle lines are read from |
| `proc as os.Process` | the decompressor, when `piped` |
| `helper as string` | FIFO directory or temp file removed on close, `""` when none |
| `piped as bool` | whether a child process feeds the source |
| `origin as string` | the path opened, or `"<stdin>"` |

### `source.Sink`

An open line sink.

| field | notes |
| --- | --- |
| `file as fs.File` | the handle lines are written to |
| `proc as os.Process` | the compressor, when `piped` |
| `helper as string` | FIFO directory removed on close, `""` when none |
| `piped as bool` | whether a child process drains the sink |
| `origin as string` | the path opened |

| function | returns | notes |
| --- | --- | --- |
| `source.isGzip(path)` | `bool` | first two bytes are `1f 8b`; raises when the file is missing |
| `source.open(path)` | `Source` | gzip detected from content |
| `source.openPlain(path)` | `Source` | raises when the file is missing |
| `source.openStdin()` | `Source` | via `/dev/stdin`, Unix only; raises when it is missing, `runtime` when it cannot be reopened |
| `source.openGzip(path)` | `Source` | streams `gzip -dc` through a FIFO; raises without `exec` or when `mkfifo` fails |
| `source.openGzipBuffered(path)` | `Source` | decompresses the whole file in memory to a temp file; small inputs only |
| `source.hasNext(s)` | `bool` | false at end of stream |
| `source.nextLine(s)` | `string` | next line, terminator stripped; raises at end of stream |
| `source.close(s)` | `int` | decompressor exit code, 0 without one; releases child and temp files even when it raises |
| `source.closeStrict(s)` | - | raises when the decompressor exited non-zero; use after reading to the end |
| `source.create(path)` | `Sink` | gzip when the name ends in `.gz` |
| `source.createPlain(path)` | `Sink` | creates or truncates |
| `source.createGzip(path)` | `Sink` | `gzip -c` reading a FIFO; raises without `exec` or when `mkfifo` fails |
| `source.writeText(k, text)` | - | writes `text` as is |
| `source.writeLine(k, line)` | - | writes `line` and `"\n"` |
| `source.closeSink(k)` | `int` | flushes, waits for the compressor; its exit code, 0 without one |

## `fastq`

FASTQ records, reader and writer ([Formats](formats.md#fastq)). Errors are
`ngs.fastq`.

| constant | value | notes |
| --- | --- | --- |
| `fastq.PHRED33` | `33` | Sanger / Illumina 1.8+ quality offset |
| `fastq.PHRED64` | `64` | legacy Illumina 1.3-1.7 offset |

### `fastq.Read`

| field | notes |
| --- | --- |
| `id as string` | header up to the first space, without `@` |
| `desc as string` | header after the first space, `""` when none |
| `seq as string` | sequence |
| `qual as string` | raw ASCII qualities, same length as `seq` |

### `fastq.Reader`

| field | notes |
| --- | --- |
| `src as source.Source` | the underlying line source |

| function | returns | notes |
| --- | --- | --- |
| `fastq.open(path)` | `Reader` | gzip transparent |
| `fastq.openStdin()` | `Reader` | standard input |
| `fastq.hasNext(r)` | `bool` | another record available |
| `fastq.next(r)` | `Read` | raises on a truncated or malformed record, naming its header |
| `fastq.close(r)` | `int` | decompressor exit code |
| `fastq.closeStrict(r)` | - | raises when the decompressor failed |
| `fastq.parse(header, seq, plus, qual)` | `Read` | raises on a bad `@` or `+` marker or a length mismatch |
| `fastq.format(rec)` | `string` | four lines, newline-terminated, separator `+` |
| `fastq.toFasta(rec)` | `string` | two lines, newline-terminated |
| `fastq.headerOf(rec)` | `string` | `id`, plus `" " + desc` when present |
| `fastq.write(k, rec)` | - | writes `format(rec)` to a sink |
| `fastq.length(rec)` | `int` | bases |
| `fastq.qualities(rec, offset)` | `list of int` | one Phred score per base; one step per base |
| `fastq.meanQuality(rec, offset)` | `float` | arithmetic mean Phred, 0.0 for an empty read |
| `fastq.minQuality(rec, offset)` | `int` | lowest Phred, -1 for an empty read |
| `fastq.gcPercent(rec)` | `float` | 0-100, 0.0 for an empty read |
| `fastq.ambiguousCount(rec)` | `int` | N bases |
| `fastq.subseq(rec, start, end)` | `Read` | 0-based half-open slice of `seq` and `qual`; raises when out of range |
| `fastq.reverseComplement(rec)` | `Read` | qualities reversed to match |

## `bed`

BED as `intervals.Interval` ([Formats](formats.md#bed)). Errors are `ngs.bed`.

| constant | value | notes |
| --- | --- | --- |
| `bed.MISSING` | `"."` | written for an absent name or strand |

### `bed.Reader`

| field | notes |
| --- | --- |
| `src as source.Source` | the underlying line source |

### `bed.Next`

| variant | notes |
| --- | --- |
| `Feature` | the next record, in `value` (`intervals.Interval`) |
| `Done` | end of stream |

| function | returns | notes |
| --- | --- | --- |
| `bed.open(path)` | `Reader` | gzip transparent |
| `bed.openStdin()` | `Reader` | standard input |
| `bed.next(r)` | `Next` | skips blank, `#`, `track` and `browser` lines |
| `bed.close(r)` | `int` | decompressor exit code |
| `bed.closeStrict(r)` | - | raises when the decompressor failed |
| `bed.parse(line)` | `intervals.Interval` | first six columns; raises on fewer than three or a non-numeric field |
| `bed.format(iv)` | `string` | BED6, no newline; `.` for an empty name or strand |
| `bed.readAll(path)` | `list of intervals.Interval` | whole file, file order |
| `bed.writeAll(path, ivs)` | - | gzip when the name ends in `.gz`; raises when the compressor fails |
| `bed.write(k, iv)` | - | one BED6 line to an open sink |

## `gff`

GFF3 and GTF ([Formats](formats.md#gff3-and-gtf)). Errors are `ngs.gff`.

| constant | value | notes |
| --- | --- | --- |
| `gff.MISSING` | `"."` | placeholder for an absent column |
| `gff.STYLE_GFF3` | `"gff3"` | `key=value;` attributes, percent-escaped |
| `gff.STYLE_GTF` | `"gtf"` | `key "value";` attributes |

### `gff.Feature`

| field | notes |
| --- | --- |
| `seqid as string` | reference sequence |
| `origin as string` | the `source` column |
| `kind as string` | the `type` column |
| `start as int` | 1-based, inclusive, as written |
| `end as int` | 1-based, inclusive, as written |
| `score as float` | 0.0 when absent |
| `strand as string` | `+`, `-` or `.` |
| `phase as int` | CDS phase 0-2, -1 when absent |
| `attributes as map of string to string` | the ninth column, in file order |

### `gff.Reader`

| field | notes |
| --- | --- |
| `src as source.Source` | the underlying line source |

### `gff.Next`

| variant | notes |
| --- | --- |
| `Record` | the next feature, in `value` (`Feature`) |
| `Done` | end of stream, or a `##FASTA` line |

| function | returns | notes |
| --- | --- | --- |
| `gff.open(path)` | `Reader` | gzip transparent |
| `gff.openStdin()` | `Reader` | standard input |
| `gff.next(r)` | `Next` | skips blank and `#` lines, stops at `##FASTA` |
| `gff.close(r)` | `int` | decompressor exit code |
| `gff.closeStrict(r)` | - | raises when the decompressor failed |
| `gff.parse(line)` | `Feature` | style detected per record; raises on fewer than nine columns |
| `gff.parseAttributes(text)` | `map of string to string` | GFF3 values percent-decoded; empty for `""` or `.` |
| `gff.styleOf(text)` | `string` | `STYLE_GFF3` or `STYLE_GTF`, from the first pair |
| `gff.attribute(f, key, fallback)` | `string` | value, or `fallback` |
| `gff.identifier(f)` | `string` | `ID`, else `transcript_id`, else `gene_id`, else `""` |
| `gff.name(f)` | `string` | `Name`, else `gene_name`, else `identifier(f)` |
| `gff.toInterval(f)` | `intervals.Interval` | `start - 1`, `end` kept, named by `name(f)` |
| `gff.fromInterval(iv, origin, kind)` | `Feature` | `start + 1`; `ID` from the interval's name; phase -1 |
| `gff.format(f, style)` | `string` | nine columns, no newline; score 0.0 written as `.` |
| `gff.formatAttributes(attrs, style)` | `string` | `.` for an empty map; GFF3 escapes only reserved characters |
| `gff.readAll(path)` | `list of Feature` | whole file up to `##FASTA`; closes without checking the decompressor |
| `gff.readIntervals(path, kind)` | `list of intervals.Interval` | converted features of type `kind`, `""` for all |

## `sam`

SAM alignments ([Formats](formats.md#sam)). Errors are `ngs.sam`.

| constant | value | notes |
| --- | --- | --- |
| `sam.FLAG_PAIRED` | `0x1` | paired in sequencing |
| `sam.FLAG_PROPER_PAIR` | `0x2` | each segment properly aligned |
| `sam.FLAG_UNMAPPED` | `0x4` | segment unmapped |
| `sam.FLAG_MATE_UNMAPPED` | `0x8` | mate unmapped |
| `sam.FLAG_REVERSE` | `0x10` | segment on the reverse strand |
| `sam.FLAG_MATE_REVERSE` | `0x20` | mate on the reverse strand |
| `sam.FLAG_FIRST` | `0x40` | first segment of the template |
| `sam.FLAG_LAST` | `0x80` | last segment of the template |
| `sam.FLAG_SECONDARY` | `0x100` | secondary alignment |
| `sam.FLAG_QC_FAIL` | `0x200` | failed platform or vendor QC |
| `sam.FLAG_DUPLICATE` | `0x400` | PCR or optical duplicate |
| `sam.FLAG_SUPPLEMENTARY` | `0x800` | supplementary alignment |
| `sam.MISSING` | `"*"` | placeholder for an absent column |

### `sam.Record`

| field | notes |
| --- | --- |
| `qname as string` | query template name |
| `flag as int` | FLAG bits |
| `rname as string` | reference name, `*` when unmapped |
| `pos as int` | 1-based leftmost position, 0 when unmapped |
| `mapq as int` | mapping quality |
| `cigar as string` | raw CIGAR, `*` when unavailable |
| `rnext as string` | mate's reference name |
| `pnext as int` | mate's 1-based position |
| `tlen as int` | observed template length, signed |
| `seq as string` | segment sequence, `*` when absent |
| `qual as string` | raw ASCII qualities, `*` when absent |
| `tags as list of string` | optional `TAG:TYPE:VALUE` columns, verbatim |

### `sam.CigarOp`

| field | notes |
| --- | --- |
| `length as int` | bases covered |
| `op as string` | one of `M I D N S H P = X` |

### `sam.Reader`

| field | notes |
| --- | --- |
| `src as source.Source` | the underlying line source |
| `header as list of string` | the `@` lines, empty on stdin |

### `sam.Next`

| variant | notes |
| --- | --- |
| `Alignment` | the next record, in `value` (`Record`) |
| `Done` | end of stream |

| function | returns | notes |
| --- | --- | --- |
| `sam.open(path)` | `Reader` | reads the header in a separate pass |
| `sam.openStdin()` | `Reader` | standard input; `header` empty |
| `sam.next(r)` | `Next` | skips `@` and blank lines |
| `sam.close(r)` | `int` | decompressor exit code |
| `sam.closeStrict(r)` | - | raises when the decompressor failed |
| `sam.parse(line)` | `Record` | raises on fewer than eleven columns or a non-numeric field |
| `sam.format(rec)` | `string` | one record, no newline |
| `sam.hasFlag(rec, bit)` | `bool` | true when any bit of `bit` is set |
| `sam.isUnmapped(rec)` | `bool` | FLAG `0x4` |
| `sam.isReverse(rec)` | `bool` | FLAG `0x10` |
| `sam.isSecondary(rec)` | `bool` | FLAG `0x100` |
| `sam.isSupplementary(rec)` | `bool` | FLAG `0x800` |
| `sam.isDuplicate(rec)` | `bool` | FLAG `0x400` |
| `sam.isPrimary(rec)` | `bool` | none of `0x4`, `0x100`, `0x800` |
| `sam.strandOf(rec)` | `string` | `+`, `-`, or `.` when unmapped |
| `sam.parseCigar(cigar)` | `list of CigarOp` | empty for `*`; raises on a malformed CIGAR |
| `sam.referenceLength(cigar)` | `int` | sum of `M D N = X` |
| `sam.queryLength(cigar)` | `int` | sum of `M I S = X` |
| `sam.toInterval(rec)` | `intervals.Interval` | `pos - 1`, end from CIGAR (sequence length for `*`); named by `qname`, scored by MAPQ; raises when unmapped |
| `sam.tagValue(rec, tagName, fallback)` | `string` | value with `TAG:TYPE:` removed, or `fallback` |
| `sam.referenceSizes(r)` | `map of string to int` | `@SQ` `SN` to `LN` |
| `sam.headerOf(path)` | `list of string` | the `@` lines only |

## `vcf`

VCF variants ([Formats](formats.md#vcf)). Errors are `ngs.vcf`.

| constant | value | notes |
| --- | --- | --- |
| `vcf.MISSING` | `"."` | placeholder for an absent column |
| `vcf.FIXED_COLUMNS` | `9` | columns before the first sample: eight fixed plus FORMAT |

### `vcf.Variant`

| field | notes |
| --- | --- |
| `chrom as string` | reference sequence |
| `pos as int` | 1-based position of the first REF base |
| `id as string` | `""` when `.` |
| `ref as string` | reference allele |
| `alt as list of string` | alternate alleles, empty when `.` |
| `qual as float` | Phred-scaled quality, -1.0 when `.` |
| `filter as string` | FILTER verbatim |
| `info as map of string to string` | INFO; a flag maps to `""` |
| `format as list of string` | FORMAT keys in order |
| `samples as list of string` | per-sample columns verbatim |

### `vcf.Reader`

| field | notes |
| --- | --- |
| `src as source.Source` | the underlying line source |
| `meta as list of string` | `##` lines, empty on stdin |
| `samples as list of string` | sample names from `#CHROM`, empty on stdin |

### `vcf.Next`

| variant | notes |
| --- | --- |
| `Record` | the next variant, in `value` (`Variant`) |
| `Done` | end of stream |

| function | returns | notes |
| --- | --- | --- |
| `vcf.open(path)` | `Reader` | reads the header in a separate pass |
| `vcf.openStdin()` | `Reader` | standard input; `meta` and `samples` empty |
| `vcf.next(r)` | `Next` | skips `#` and blank lines |
| `vcf.close(r)` | `int` | decompressor exit code |
| `vcf.closeStrict(r)` | - | raises when the decompressor failed |
| `vcf.headerOf(path)` | `list of string` | the leading `#` lines |
| `vcf.parse(line)` | `Variant` | raises on fewer than eight columns or a non-numeric POS or QUAL |
| `vcf.parseInfo(text)` | `map of string to string` | value from the first `=`; flags map to `""`; empty for `.` |
| `vcf.info(v, key, fallback)` | `string` | INFO value, or `fallback` |
| `vcf.hasInfo(v, key)` | `bool` | key present |
| `vcf.passed(v)` | `bool` | FILTER is `PASS`, `.` or empty |
| `vcf.isSnv(v)` | `bool` | REF and every ALT one of `A C G T N` |
| `vcf.isIndel(v)` | `bool` | an ALT of different length from REF; symbolic and breakend ALTs excluded |
| `vcf.sampleValue(v, sampleIndex, key, fallback)` | `string` | FORMAT value for a 0-based sample, or `fallback` |
| `vcf.genotype(v, sampleIndex)` | `string` | GT, `./.` when unavailable |
| `vcf.toInterval(v)` | `intervals.Interval` | `pos - 1`, spans REF; named by ID; QUAL as score, 0.0 when absent |
| `vcf.format(v)` | `string` | one record, no newline; QUAL rendered as a float |
| `vcf.formatInfo(entries)` | `string` | `key=value` or bare flag, `;`-joined; `.` when empty |

## `intervals`

Interval arithmetic, 0-based half-open ([Interval arithmetic](intervals.md)).
Errors are `ngs.intervals`. Set operations return results in the order of
`sort(a)`.

| constant | value | notes |
| --- | --- | --- |
| `intervals.KEY_WIDTH` | `12` | digits coordinates are zero-padded to in sort keys |
| `intervals.MAX_COORDINATE` | `1000000000000` | one past the largest allowed coordinate (10^12) |

### `intervals.Interval`

| field | notes |
| --- | --- |
| `chrom as string` | reference sequence |
| `start as int` | 0-based, inclusive |
| `end as int` | 0-based, exclusive |
| `name as string` | `""` when unnamed |
| `score as float` | 0.0 when absent |
| `strand as string` | `+`, `-` or `.` |

| function | returns | notes |
| --- | --- | --- |
| `intervals.make(chrom, start, end)` | `Interval` | unnamed, strand `.`; validated as `named` |
| `intervals.named(chrom, start, end, name, score, strand)` | `Interval` | raises on a negative or inverted range or an end of `MAX_COORDINATE` or more |
| `intervals.length(iv)` | `int` | `end - start` |
| `intervals.overlaps(a, b)` | `bool` | share a base; book-ended intervals do not |
| `intervals.overlapLength(a, b)` | `int` | shared bases, 0 when none |
| `intervals.contains(outer, inner)` | `bool` | same chromosome, `inner` within `outer` |
| `intervals.intersection(a, b)` | `Interval` | shared region with `a`'s name, score and strand; raises when they do not overlap |
| `intervals.sortKey(iv)` | `string` | chrom, start, end, tab-separated, coordinates zero-padded; exported as the `func` value for `lists.sortBy` |
| `intervals.sort(ivs)` | `list of Interval` | by chromosome as text, start, end; sorted input returned as is |
| `intervals.merge(ivs, maxGap)` | `list of Interval` | per chromosome; 0 joins overlapping and book-ended (`bedtools merge -d 0`); names, scores, strands cleared |
| `intervals.totalLength(ivs)` | `int` | bases in the union |
| `intervals.intersect(a, b)` | `list of Interval` | one per overlapping pair, `a`'s name (`bedtools intersect`) |
| `intervals.countOverlaps(a, b)` | `list of int` | targets overlapping each query (`intersect -c`) |
| `intervals.coverage(a, b)` | `list of int` | bases of each query covered by merged `b` |
| `intervals.subtract(a, b)` | `list of Interval` | uncovered pieces of each query, names kept (`bedtools subtract`) |
| `intervals.complement(ivs, sizes)` | `list of Interval` | gaps on the chromosomes in `sizes`, sorted (`bedtools complement`) |
| `intervals.slop(iv, left, right, limit)` | `Interval` | grown, clamped at 0 and at `limit` when `limit > 0` |
| `intervals.byChromosome(ivs)` | `map of string to list of Interval` | sorted groups |
| `intervals.toBed3(iv)` | `string` | `chrom\tstart\tend`, no newline |

## `qc`

One-pass read QC ([Read quality control](qc.md)).

| constant | value | notes |
| --- | --- | --- |
| `qc.MAX_PHRED` | `94` | size of the quality histogram, Phred 0-93 |
| `qc.Q30` | `30` | threshold for `q30Percent` |
| `qc.DEFAULT_BATCH` | `2000` | records per counting batch |
| `qc.DEFAULT_BATCH_BASES` | `4000000` | bases per counting batch |
| `qc.ADAPTER_TRUSEQ` | `"AGATCGGAAGAGC"` | Illumina TruSeq / NEBNext universal |
| `qc.ADAPTER_NEXTERA` | `"CTGTCTCTTATACACATCT"` | Nextera / transposase |
| `qc.ADAPTER_SMALL_RNA` | `"TGGAATTCTCGG"` | Illumina small-RNA 3' |
| `qc.ARTEFACT_POLY_G` | `"GGGGGGGGGG"` | poly-G of two-colour chemistry (NovaSeq, NextSeq) |

### `qc.Options`

The zero value is valid.

| field | notes |
| --- | --- |
| `maxReads as int` | stop after sampling this many; 0 is no limit |
| `everyNth as int` | sample one read in N; 0 or 1 is every read |
| `perBase as bool` | collect per-position quality, one step per base |
| `offset as int` | quality offset; 0 means 33 |
| `batchSize as int` | records per batch; 0 means `DEFAULT_BATCH` |
| `maxBatchBases as int` | bases per batch; 0 means `DEFAULT_BATCH_BASES` |
| `adapters as list of string` | probes to count; empty means `defaultAdapters()` |

### `qc.Report`

| field | notes |
| --- | --- |
| `reads as int` | records read |
| `sampled as int` | records measured |
| `bases as int` | bases in sampled records |
| `minLength as int` | shortest sampled read, 0 when none |
| `maxLength as int` | longest sampled read |
| `gcBases as int` | G and C bases |
| `nBases as int` | N bases |
| `qualHist as list of int` | bases per Phred score, indexed 0-93 |
| `lengthHist as map of int to int` | read length to read count |
| `adapterHits as map of string to int` | probe to occurrences |
| `perBaseSum as list of int` | quality sum per position, `perBase` only |
| `perBaseCount as list of int` | reads reaching each position, `perBase` only |

| function | returns | notes |
| --- | --- | --- |
| `qc.options()` | `Options` | every read, Phred+33, default batch bounds |
| `qc.defaultAdapters()` | `list of string` | the four probe constants |
| `qc.analyze(path, opts)` | `Report` | one streaming pass; raises `ngs.qc` on qualities outside Phred 0-93 |
| `qc.meanQuality(rep)` | `float` | mean Phred over sampled bases |
| `qc.q30Percent(rep)` | `float` | percentage of bases at Q30 or above |
| `qc.gcPercent(rep)` | `float` | percentage G+C |
| `qc.meanLength(rep)` | `float` | mean sampled read length |
| `qc.perBaseMean(rep)` | `list of float` | mean quality per position, empty without `perBase` |
| `qc.summary(rep)` | `string` | multi-line text report |

## `trim`

Trimming, filtering and subsampling ([Trimming and filtering](preprocessing.md)).
Errors are `ngs.trim`.

| constant | value | notes |
| --- | --- | --- |
| `trim.DEFAULT_WINDOW` | `4` | sliding window width |
| `trim.DEFAULT_OVERLAP` | `6` | shortest partial 3' adapter match |

### `trim.Options`

The zero value trims nothing; its `maxAmbiguous` of 0 still drops reads
holding an N.

| field | notes |
| --- | --- |
| `trimFront as int` | bases cut from the 5' end |
| `trimTail as int` | bases cut from the 3' end |
| `adapters as list of string` | adapters to clip; empty disables |
| `minOverlap as int` | shortest partial adapter; 0 means `DEFAULT_OVERLAP` |
| `qualityThreshold as int` | window mean quality to keep; 0 disables |
| `windowSize as int` | window width; 0 means `DEFAULT_WINDOW` |
| `maxLength as int` | length cap; 0 is none |
| `minLength as int` | shorter reads dropped after trimming |
| `maxAmbiguous as int` | reads with more N bases dropped; -1 disables |
| `minMeanQuality as int` | reads below this mean dropped; 0 disables |
| `offset as int` | quality offset; 0 means 33 |

### `trim.Stats`

| field | notes |
| --- | --- |
| `readsIn as int` | records read |
| `readsOut as int` | records written |
| `basesIn as int` | bases read |
| `basesOut as int` | bases written |
| `adapterTrimmed as int` | records shortened by adapter clipping |
| `qualityTrimmed as int` | records shortened by quality trimming |
| `tooShort as int` | dropped by `minLength` |
| `tooAmbiguous as int` | dropped by `maxAmbiguous` |
| `lowQuality as int` | dropped by `minMeanQuality` |

| function | returns | notes |
| --- | --- | --- |
| `trim.options()` | `Options` | built-in adapters, Q20 window of 4, `minLength` 25, `maxAmbiguous` -1 |
| `trim.run(inPath, outPath, opts)` | `Stats` | streams; raises when input and output are the same file or the compressor fails |
| `trim.trimRead(rec, opts)` | `fastq.Read` | the four trimming steps on one record |
| `trim.clipFixed(rec, front, tail)` | `fastq.Read` | empty when the cuts overlap |
| `trim.clipAdapters(rec, adapters, minOverlap)` | `fastq.Read` | cut at the earliest full or 3' partial adapter |
| `trim.adapterStart(seq, adapters, minOverlap)` | `int` | 0-based cut position, -1 when none; `minOverlap` below 1 means 6 |
| `trim.clipQuality(rec, threshold, window, offset)` | `fastq.Read` | cut at the first window below `threshold`; 0 disables |
| `trim.qualityKeepLength(qual, threshold, window, offset)` | `int` | bases kept from the 5' end |
| `trim.truncate(rec, maxLength)` | `fastq.Read` | 5' end kept; 0 or less is no cap |
| `trim.keeps(rec, opts)` | `bool` | passes length, N and mean-quality filters |
| `trim.subsampleEvery(inPath, outPath, everyNth, maxReads)` | `Stats` | every Nth record; `maxReads` 0 is no limit |
| `trim.subsampleFraction(inPath, outPath, fraction)` | `Stats` | random, `fraction` in (0, 1]; seed with `math.randSeed` |

## `pipeline`

External tools ([Pipeline orchestration](pipelines.md)). Errors are
`ngs.pipeline`.

### `pipeline.Job`

| field | notes |
| --- | --- |
| `name as string` | label in logs and results |
| `argv as list of string` | command and arguments, run without a shell |
| `inputs as list of string` | files read |
| `outputs as list of string` | files produced |
| `stdin as string` | text fed to standard input |
| `allowFailure as bool` | a non-zero exit is reported, not raised |

### `pipeline.Result`

| field | notes |
| --- | --- |
| `name as string` | the job's label |
| `exitCode as int` | process exit code, 0 when skipped |
| `stdout as string` | captured standard output |
| `stderr as string` | captured standard error |
| `skipped as bool` | the job was up to date |
| `allowFailure as bool` | copied from the job |

| function | returns | notes |
| --- | --- | --- |
| `pipeline.job(name, argv)` | `Job` | raises on an empty argv |
| `pipeline.shell(name, cmdline)` | `Job` | runs `sh -c cmdline`; quote untrusted values with `core.shellQuote` |
| `pipeline.needing(j, inputs)` | `Job` | copy with `inputs` set |
| `pipeline.producing(j, outputs)` | `Job` | copy with `outputs` set |
| `pipeline.feeding(j, text)` | `Job` | copy with `stdin` set |
| `pipeline.tolerant(j)` | `Job` | copy with `allowFailure` set |
| `pipeline.needsRun(j)` | `bool` | true without outputs, with an output missing, or with an input not older than the oldest output; raises on a missing input |
| `pipeline.run(j)` | `Result` | skips an up-to-date job; raises without `exec` |
| `pipeline.runAll(jobs)` | `list of Result` | in order; raises at the first real failure |
| `pipeline.runParallel(jobs, width)` | `list of Result` | waves of `width` (0 or less: `os.NCPU`), input order; does not check |
| `pipeline.check(res)` | - | raises on a non-zero exit without `allowFailure` |
| `pipeline.checkAll(results)` | - | `check` on each, raising at the first failure |
| `pipeline.allOk(results)` | `bool` | no real failure |
| `pipeline.hasTool(name)` | `bool` | `core.hasTool` |
| `pipeline.requireTools(names)` | - | raises listing every missing tool |
| `pipeline.workDir(tag)` | `string` | fresh temp directory; the caller removes it |
| `pipeline.cleanUp(dir)` | - | recursive delete; raises for a root, the working directory, `.` or `..` |
| `pipeline.summary(results)` | `string` | one `skip` / `ok` / `warn` / `fail` line per job |

## `ngs`

The entry module. Errors are `ngs.ngs`.

| constant | value | notes |
| --- | --- | --- |
| `ngs.VERSION` | `"0.1.0"` | the deck's version |
| `ngs.NAME` | `"@deckyard/ngs"` | the package name |

### `ngs.Capabilities`

| field | notes |
| --- | --- |
| `streamingGzip as bool` | the host has `exec`, which streaming gzip needs (`gzipOnPath` is not consulted); false means `source.openGzipBuffered` |
| `externalTools as bool` | `pipeline` can run programs |
| `gzipOnPath as bool` | a `gzip` executable is on `PATH` |

| function | returns | notes |
| --- | --- | --- |
| `ngs.capabilities()` | `Capabilities` | what this interpreter and host support |
| `ngs.describe()` | `string` | name, version and capabilities on one line |
| `ngs.requireFullHost()` | - | raises without the `exec` capability |
