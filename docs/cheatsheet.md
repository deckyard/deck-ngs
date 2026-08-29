# Cheatsheet

Every exported name in the deck, by module. A function's return type is not
declared in Jennifer, so it is given in the description.

For the reasoning behind these calls, start from the topic pages - this page is
the index, not the explanation.


## `core`

Shared primitives the rest of the deck is built on: error construction, Go-speed counting, tab-delimited field access.


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `core.ERROR_KIND` | `string` | `"ngs"` | Deck-wide error kind prefix, so a caller can match on `$e.kind`. |
| `core.TAB` | `string` | `"\t"` | Tab, the field separator every text format in this deck uses. |


### Functions

| function | description |
| --- | --- |
| `core.fail(what as string, message as string)` | Raise a deck error with a stable kind. |
| `core.countIn(haystack as string, needle as string)` | Count non-overlapping occurrences of `needle` in `haystack`. |
| `core.parseInt(what as string, label as string, raw as string)` | Parse an integer field, reporting a malformed value as a deck error. |
| `core.parseFloat(what as string, label as string, raw as string)` | Parse a float field, reporting a malformed value as a deck error. |
| `core.reraise(what as string, line as string, e as Error)` | Re-raise a parse failure as a deck error naming the record it came from. |
| `core.splitTab(line as string)` | Split a tab-delimited line into its fields. |
| `core.field(fields as list of string, i as int, fallback as string)` | Field `i` of `fields`, or `fallback` when the column is absent. |
| `core.intField(what as string, fields as list of string, i as int, fallback as int)` | Integer field `i`, or `fallback` when absent, empty or a placeholder dot. |
| `core.floatField(what as string, fields as list of string, i as int, fallback as float)` | Float field `i`, or `fallback` when absent, empty or a placeholder dot. |
| `core.isSkippable(line as string, marker as string)` | Whether a line should be skipped by a record reader: empty or a comment. |
| `core.padNum(n as int, width as int)` | Left-pad an integer with zeros to `width`, for building sortable string keys. |
| `core.shellQuote(s as string)` | Quote a string for safe inclusion in a POSIX `sh -c` command line. |
| `core.reverseComplement(seq as string)` | Reverse complement of a DNA sequence. |
| `core.complementBase(b as string)` | Complement of a single nucleotide character. |
| `core.gcCount(seq as string)` | G+C count of a sequence, counted at Go speed. |
| `core.nCount(seq as string)` | Ambiguous-base (N) count of a sequence, counted at Go speed. |
| `core.hasTool(name as string)` | Whether an executable is on `PATH`. |
| `core.sameFile(a as string, b as string)` | Whether two paths resolve to the same file. |


## `source`

Streaming line sources and sinks, with transparent gzip in both directions.


### Types

| type | kind | description |
| --- | --- | --- |
| `source.Source` | struct | An open line source. Treat it as an opaque handle and pass it around by value: the `fs.File` and `os.Process` fields are registry handles, so copies share the same underlying stream. |
| `source.Sink` | struct | An open line sink, the write-side counterpart of `Source`. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `source.GZIP_MAGIC1` | `int` | `0x1f` | First byte of the gzip magic number (RFC 1952). |
| `source.GZIP_MAGIC2` | `int` | `0x8b` | Second byte of the gzip magic number (RFC 1952). |


### Functions

| function | description |
| --- | --- |
| `source.isGzip(path as string)` | Whether a file starts with the gzip magic number. |
| `source.open(path as string)` | Open a file for streaming line reads, decompressing gzip transparently. |
| `source.openPlain(path as string)` | Open an uncompressed file for streaming line reads. |
| `source.openStdin()` | Open standard input as a line source, so the deck composes with Unix pipes. |
| `source.openGzip(path as string)` | Open a gzipped file for streaming line reads at constant memory. |
| `source.openGzipBuffered(path as string)` | Open a gzipped file by decompressing it whole, for hosts without `exec`. |
| `source.hasNext(s as Source)` | Whether another line can be read from the source. |
| `source.nextLine(s as Source)` | Read the next line, with its line terminator stripped. |
| `source.close(s as Source)` | Close a source and release its child process and temporary files. |
| `source.closeStrict(s as Source)` | Close a fully-consumed source, raising if its decompressor failed. |
| `source.create(path as string)` | Open a file for streaming line writes, gzip-compressing when the name ends in `.gz`. |
| `source.createPlain(path as string)` | Open an uncompressed file for streaming line writes. |
| `source.createGzip(path as string)` | Open a gzip-compressed file for streaming line writes at constant memory. |
| `source.writeText(k as Sink, text as string)` | Write text to a sink verbatim, adding no line terminator. |
| `source.writeLine(k as Sink, line as string)` | Write one line to a sink, appending a newline. |
| `source.closeSink(k as Sink)` | Close a sink, flush it and wait for its compressor to finish. |


## `fastq`

FASTQ records, streaming reader and writer, per-record accessors.


### Types

| type | kind | description |
| --- | --- | --- |
| `fastq.Read` | struct | One FASTQ record. |
| `fastq.Reader` | struct | A streaming FASTQ cursor. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `fastq.PHRED33` | `int` | `33` | Sanger / Illumina 1.8+ quality offset; the modern default. |
| `fastq.PHRED64` | `int` | `64` | Legacy Illumina 1.3-1.7 quality offset. |


### Functions

| function | description |
| --- | --- |
| `fastq.open(path as string)` | Open a FASTQ file for streaming, decompressing gzip transparently. |
| `fastq.openStdin()` | Open standard input as a FASTQ stream, for use in Unix pipes. |
| `fastq.hasNext(r as Reader)` | Whether another record is available. |
| `fastq.next(r as Reader)` | Read the next record. |
| `fastq.close(r as Reader)` | Close the reader and release its decompressor. |
| `fastq.closeStrict(r as Reader)` | Close a fully-consumed reader, raising if its decompressor failed. |
| `fastq.parse(header as string, seq as string, plus as string, qual as string)` | Build a record from its four raw lines, validating the FASTQ shape. |
| `fastq.format(rec as Read)` | Render a record back to its four-line FASTQ text, newline-terminated. |
| `fastq.toFasta(rec as Read)` | Render a record as two-line FASTA text, newline-terminated. |
| `fastq.headerOf(rec as Read)` | The full header text of a record: id plus description when present. |
| `fastq.write(k as source.Sink, rec as Read)` | Write a record to a sink as FASTQ. |
| `fastq.length(rec as Read)` | Read length in bases. |
| `fastq.qualities(rec as Read, offset as int)` | Decode the quality string to Phred scores. |
| `fastq.meanQuality(rec as Read, offset as int)` | Mean Phred quality of a record. |
| `fastq.minQuality(rec as Read, offset as int)` | Lowest Phred score in a record. |
| `fastq.gcPercent(rec as Read)` | G+C content of a record as a percentage, counted at Go speed. |
| `fastq.ambiguousCount(rec as Read)` | Number of ambiguous (N) bases in a record, counted at Go speed. |
| `fastq.subseq(rec as Read, start as int, end as int)` | A sub-range of a record, slicing sequence and quality together. |
| `fastq.reverseComplement(rec as Read)` | Reverse complement of a record, with the quality string reversed to match. |


## `bed`

Streaming BED reading and writing, straight onto the interval engine.


### Types

| type | kind | description |
| --- | --- | --- |
| `bed.Reader` | struct | A streaming BED cursor. |
| `bed.Next` | enum | The result of advancing a BED cursor: either a feature or end-of-stream. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `bed.MISSING` | `string` | `"."` | Placeholder written for an absent name or strand. |


### Functions

| function | description |
| --- | --- |
| `bed.open(path as string)` | Open a BED file for streaming, decompressing gzip transparently. |
| `bed.openStdin()` | Open standard input as a BED stream. |
| `bed.next(r as Reader)` | Advance the cursor to the next feature, skipping blank, comment, `track` and `browser` lines. |
| `bed.close(r as Reader)` | Close the cursor and release its decompressor. |
| `bed.closeStrict(r as Reader)` | Close a fully-consumed cursor, raising if its decompressor failed. |
| `bed.parse(line as string)` | Parse one BED line into an interval. |
| `bed.format(iv as intervals.Interval)` | Render an interval as a six-column BED line. |
| `bed.readAll(path as string)` | Read a whole BED file into memory. |
| `bed.writeAll(path as string, ivs as list of intervals.Interval)` | Write intervals to a BED file, gzipping when the name ends in `.gz`. |
| `bed.write(k as source.Sink, iv as intervals.Interval)` | Write one interval to an already-open sink. |


## `gff`

Streaming GFF3 and GTF reading and writing, with attribute parsing in both styles.


### Types

| type | kind | description |
| --- | --- | --- |
| `gff.Feature` | struct | One GFF3 or GTF feature. |
| `gff.Reader` | struct | A streaming GFF/GTF cursor. |
| `gff.Next` | enum | The result of advancing a GFF cursor. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `gff.MISSING` | `string` | `"."` | Placeholder used for an absent column. |
| `gff.STYLE_GFF3` | `string` | `"gff3"` | Attribute style: `key=value;` pairs, percent-escaped. |
| `gff.STYLE_GTF` | `string` | `"gtf"` | Attribute style: `key "value";` pairs. |


### Functions

| function | description |
| --- | --- |
| `gff.open(path as string)` | Open a GFF3 or GTF file for streaming, decompressing gzip transparently. |
| `gff.openStdin()` | Open standard input as a GFF/GTF stream. |
| `gff.next(r as Reader)` | Advance the cursor to the next feature. |
| `gff.close(r as Reader)` | Close the cursor and release its decompressor. |
| `gff.closeStrict(r as Reader)` | Close a fully-consumed cursor, raising if its decompressor failed. |
| `gff.parse(line as string)` | Parse one GFF3 or GTF record. |
| `gff.parseAttributes(text as string)` | Parse a ninth-column attribute string in either GFF3 or GTF style. |
| `gff.styleOf(text as string)` | Which attribute style a ninth column uses. |
| `gff.attribute(f as Feature, key as string, fallback as string)` | Read an attribute, or a fallback when it is absent. |
| `gff.identifier(f as Feature)` | The feature's identifier: GFF3 `ID`, else GTF `transcript_id`, else `gene_id`. |
| `gff.name(f as Feature)` | The feature's display name: GFF3 `Name`, else `gene_name`, else `identifier`. |
| `gff.toInterval(f as Feature)` | Convert a feature to a half-open interval for the interval engine. |
| `gff.fromInterval(iv as intervals.Interval, origin as string, kind as string)` | Build a feature from a half-open interval, converting back to 1-based. |
| `gff.format(f as Feature, style as string)` | Render a feature as a nine-column record. |
| `gff.formatAttributes(attrs as map of string to string, style as string)` | Render an attribute map in GFF3 or GTF style. |
| `gff.readAll(path as string)` | Read a whole annotation file into memory. |
| `gff.readIntervals(path as string, kind as string)` | Read an annotation file directly as half-open intervals. |


## `sam`

Streaming SAM reading and writing, FLAG and CIGAR handling.


### Types

| type | kind | description |
| --- | --- | --- |
| `sam.Record` | struct | One SAM alignment record. |
| `sam.CigarOp` | struct | One CIGAR operation. |
| `sam.Reader` | struct | A streaming SAM cursor, carrying the file's header block. |
| `sam.Next` | enum | The result of advancing a SAM cursor. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `sam.FLAG_PAIRED` | `int` | `0x1` | Read is paired in sequencing. |
| `sam.FLAG_PROPER_PAIR` | `int` | `0x2` | Each segment is properly aligned according to the aligner. |
| `sam.FLAG_UNMAPPED` | `int` | `0x4` | This segment is unmapped. |
| `sam.FLAG_MATE_UNMAPPED` | `int` | `0x8` | The mate segment is unmapped. |
| `sam.FLAG_REVERSE` | `int` | `0x10` | This segment is on the reverse strand. |
| `sam.FLAG_MATE_REVERSE` | `int` | `0x20` | The mate segment is on the reverse strand. |
| `sam.FLAG_FIRST` | `int` | `0x40` | This is the first segment of the template. |
| `sam.FLAG_LAST` | `int` | `0x80` | This is the last segment of the template. |
| `sam.FLAG_SECONDARY` | `int` | `0x100` | This is a secondary alignment. |
| `sam.FLAG_QC_FAIL` | `int` | `0x200` | The read failed platform or vendor quality checks. |
| `sam.FLAG_DUPLICATE` | `int` | `0x400` | The read is a PCR or optical duplicate. |
| `sam.FLAG_SUPPLEMENTARY` | `int` | `0x800` | This is a supplementary alignment. |
| `sam.MISSING` | `string` | `"*"` | Placeholder used for an absent column. |


### Functions

| function | description |
| --- | --- |
| `sam.open(path as string)` | Open a SAM file for streaming, capturing its header block. |
| `sam.openStdin()` | Open standard input as a SAM stream. |
| `sam.next(r as Reader)` | Advance the cursor to the next alignment record. |
| `sam.close(r as Reader)` | Close the cursor and release its decompressor. |
| `sam.closeStrict(r as Reader)` | Close a fully-consumed cursor, raising if its decompressor failed. |
| `sam.parse(line as string)` | Parse one SAM alignment line. |
| `sam.format(rec as Record)` | Render an alignment back to a SAM record line. |
| `sam.hasFlag(rec as Record, bit as int)` | Whether a FLAG bit is set on a record. |
| `sam.isUnmapped(rec as Record)` | Whether the record is unmapped. |
| `sam.isReverse(rec as Record)` | Whether the record is mapped to the reverse strand. |
| `sam.isSecondary(rec as Record)` | Whether the record is a secondary alignment. |
| `sam.isSupplementary(rec as Record)` | Whether the record is a supplementary alignment. |
| `sam.isDuplicate(rec as Record)` | Whether the record is flagged as a PCR or optical duplicate. |
| `sam.isPrimary(rec as Record)` | Whether the record is a primary alignment: mapped, not secondary, not supplementary. |
| `sam.strandOf(rec as Record)` | The strand a record aligns to. |
| `sam.parseCigar(cigar as string)` | Split a CIGAR string into its operations. |
| `sam.referenceLength(cigar as string)` | Number of reference bases a CIGAR string consumes. |
| `sam.queryLength(cigar as string)` | Number of query bases a CIGAR string consumes. |
| `sam.toInterval(rec as Record)` | The reference interval an alignment covers, as a half-open interval. |
| `sam.tagValue(rec as Record, tagName as string, fallback as string)` | Read the value of an optional tag. |
| `sam.referenceSizes(r as Reader)` | The reference sequences declared in the header, as chromosome sizes. |
| `sam.headerOf(path as string)` | Read just the `@` header block of a SAM file. |


## `vcf`

Streaming VCF reading and writing, INFO and per-sample FORMAT access.


### Types

| type | kind | description |
| --- | --- | --- |
| `vcf.Variant` | struct | One VCF variant record. |
| `vcf.Reader` | struct | A streaming VCF cursor, carrying the file's header. |
| `vcf.Next` | enum | The result of advancing a VCF cursor. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `vcf.MISSING` | `string` | `"."` | Placeholder used for an absent column. |
| `vcf.FIXED_COLUMNS` | `int` | `9` | Number of fixed columns before the per-sample genotype columns. |


### Functions

| function | description |
| --- | --- |
| `vcf.open(path as string)` | Open a VCF file for streaming, capturing its header. |
| `vcf.openStdin()` | Open standard input as a VCF stream. |
| `vcf.next(r as Reader)` | Advance the cursor to the next variant record. |
| `vcf.close(r as Reader)` | Close the cursor and release its decompressor. |
| `vcf.closeStrict(r as Reader)` | Close a fully-consumed cursor, raising if its decompressor failed. |
| `vcf.headerOf(path as string)` | Read just the header lines of a VCF file. |
| `vcf.parse(line as string)` | Parse one VCF variant line. |
| `vcf.parseInfo(text as string)` | Parse a VCF INFO column into a map. |
| `vcf.info(v as Variant, key as string, fallback as string)` | Read an INFO field, or a fallback when it is absent. |
| `vcf.hasInfo(v as Variant, key as string)` | Whether an INFO flag is present. |
| `vcf.passed(v as Variant)` | Whether the variant passed every filter. |
| `vcf.isSnv(v as Variant)` | Whether the variant is a single-nucleotide substitution. |
| `vcf.isIndel(v as Variant)` | Whether the variant is an insertion or deletion. |
| `vcf.sampleValue(v as Variant, sampleIndex as int, key as string, fallback as string)` | The per-sample value of one FORMAT key. |
| `vcf.genotype(v as Variant, sampleIndex as int)` | The genotype call of one sample. |
| `vcf.toInterval(v as Variant)` | The reference span of a variant, as a half-open interval. |
| `vcf.format(v as Variant)` | Render a variant back to a VCF record line. |
| `vcf.formatInfo(entries as map of string to string)` | Render an INFO map back to its column text. |


## `intervals`

Interval arithmetic over genomic features - the bedtools-lite layer.


### Types

| type | kind | description |
| --- | --- | --- |
| `intervals.Interval` | struct | A half-open genomic interval. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `intervals.KEY_WIDTH` | `int` | `12` | Width used when zero-padding coordinates into sort keys. |
| `intervals.MAX_COORDINATE` | `int` | `1000000000000` | One past the largest coordinate an interval may carry. |


### Functions

| function | description |
| --- | --- |
| `intervals.make(chrom as string, start as int, end as int)` | Build an unnamed, unstranded interval. |
| `intervals.named(chrom as string,
 start as int,
 end as int,
 name as string,
 score as float,
 strand as string)` | Build a fully specified interval. |
| `intervals.length(iv as Interval)` | Width of an interval in bases. |
| `intervals.overlaps(a as Interval, b as Interval)` | Whether two intervals share at least one base. |
| `intervals.overlapLength(a as Interval, b as Interval)` | Number of bases two intervals share. |
| `intervals.contains(outer as Interval, inner as Interval)` | Whether `outer` fully contains `inner`. |
| `intervals.intersection(a as Interval, b as Interval)` | The shared region of two overlapping intervals. |
| `intervals.sortKey(iv as Interval)` | Sort key for one interval: chromosome, then zero-padded start and end. |
| `intervals.sort(ivs as list of Interval)` | Sort intervals by chromosome then start then end. |
| `intervals.merge(ivs as list of Interval, maxGap as int)` | Merge overlapping and nearby intervals, per chromosome. |
| `intervals.totalLength(ivs as list of Interval)` | Total number of bases covered, counting overlaps once. |
| `intervals.intersect(a as list of Interval, b as list of Interval)` | Pairwise intersection regions between two interval sets. |
| `intervals.countOverlaps(a as list of Interval, b as list of Interval)` | Count, for each query interval, how many target intervals it overlaps. |
| `intervals.coverage(a as list of Interval, b as list of Interval)` | Bases of each query interval covered by the target set. |
| `intervals.subtract(a as list of Interval, b as list of Interval)` | Remove from each query interval every base covered by the target set. |
| `intervals.complement(ivs as list of Interval, sizes as map of string to int)` | The gaps between intervals, across whole chromosomes. |
| `intervals.slop(iv as Interval, left as int, right as int, limit as int)` | Grow an interval on both sides, clamped at zero and an optional chromosome end. |
| `intervals.byChromosome(ivs as list of Interval)` | Group intervals by chromosome. |
| `intervals.toBed3(iv as Interval)` | Render an interval as a three-column BED line. |


## `qc`

One-pass streaming read QC - the FastQC-lite layer.


### Types

| type | kind | description |
| --- | --- | --- |
| `qc.Options` | struct | How to run an analysis. |
| `qc.Report` | struct | Everything one pass measured. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `qc.MAX_PHRED` | `int` | `94` | Highest Phred score tracked; covers the whole printable ASCII range. |
| `qc.Q30` | `int` | `30` | Quality score at or above which a base is counted towards Q30. |
| `qc.DEFAULT_BATCH` | `int` | `2000` | Records buffered before each batch of Go-speed counting. |
| `qc.DEFAULT_BATCH_BASES` | `int` | `4000000` | Bases buffered before a batch is counted, whichever limit is reached first. |
| `qc.ADAPTER_TRUSEQ` | `string` | `"AGATCGGAAGAGC"` | Illumina TruSeq / NEBNext universal adapter prefix. |
| `qc.ADAPTER_NEXTERA` | `string` | `"CTGTCTCTTATACACATCT"` | Nextera / Illumina transposase adapter prefix. |
| `qc.ADAPTER_SMALL_RNA` | `string` | `"TGGAATTCTCGG"` | Illumina small-RNA 3' adapter prefix. |
| `qc.ARTEFACT_POLY_G` | `string` | `"GGGGGGGGGG"` | Poly-G run, the two-colour-chemistry artefact seen on NovaSeq / NextSeq. |


### Functions

| function | description |
| --- | --- |
| `qc.options()` | Default analysis options: every read, no per-position quality, Phred+33. |
| `qc.defaultAdapters()` | The adapters and artefacts counted when `Options.adapters` is empty. |
| `qc.analyze(path as string, opts as Options)` | Run a single streaming pass over a FASTQ file. |
| `qc.meanQuality(rep as Report)` | Mean Phred quality across every sampled base. |
| `qc.q30Percent(rep as Report)` | Percentage of sampled bases at Q30 or better. |
| `qc.gcPercent(rep as Report)` | G+C content of the sampled bases as a percentage. |
| `qc.meanLength(rep as Report)` | Mean sampled read length. |
| `qc.perBaseMean(rep as Report)` | Mean quality at each read position. |
| `qc.summary(rep as Report)` | Render a report as a human-readable text summary. |


## `trim`

Trimming, filtering and subsampling - the fastp-lite layer.


### Types

| type | kind | description |
| --- | --- | --- |
| `trim.Options` | struct | How to trim and filter. |
| `trim.Stats` | struct | What a run did. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `trim.DEFAULT_WINDOW` | `int` | `4` | Default sliding window width for quality trimming. |
| `trim.DEFAULT_OVERLAP` | `int` | `6` | Shortest partial adapter match accepted at a read's 3' end. |


### Functions

| function | description |
| --- | --- |
| `trim.options()` | Default options: adapter clipping on, sliding-window quality trimming at Q20, and reads shorter than 25 bases discarded. |
| `trim.run(inPath as string, outPath as string, opts as Options)` | Stream a FASTQ file through trimming and filtering into a new file. |
| `trim.trimRead(rec as fastq.Read, opts as Options)` | Apply every configured trimming step to one record. |
| `trim.clipFixed(rec as fastq.Read, front as int, tail as int)` | Remove a fixed number of bases from each end. |
| `trim.clipAdapters(rec as fastq.Read, adapters as list of string, minOverlap as int)` | Clip the read at the first adapter occurrence, including a partial adapter running off the 3' end. |
| `trim.adapterStart(seq as string, adapters as list of string, minOverlap as int)` | The earliest position at which an adapter starts in a sequence. |
| `trim.clipQuality(rec as fastq.Read, threshold as int, window as int, offset as int)` | Trim the 3' end where a sliding window's mean quality falls below a threshold. |
| `trim.qualityKeepLength(qual as string, threshold as int, window as int, offset as int)` | How many bases survive sliding-window quality trimming. |
| `trim.truncate(rec as fastq.Read, maxLength as int)` | Cap a record's length, keeping the 5' end. |
| `trim.keeps(rec as fastq.Read, opts as Options)` | Whether a trimmed record passes the configured filters. |
| `trim.subsampleEvery(inPath as string, outPath as string, everyNth as int, maxReads as int)` | Copy every Nth record of a FASTQ file to a new file. |
| `trim.subsampleFraction(inPath as string, outPath as string, fraction as float)` | Copy a random fraction of a FASTQ file to a new file. |


## `pipeline`

Run, skip, check and fan out external tools.


### Types

| type | kind | description |
| --- | --- | --- |
| `pipeline.Job` | struct | One step of a pipeline. |
| `pipeline.Result` | struct | What running a job produced. |


### Functions

| function | description |
| --- | --- |
| `pipeline.job(name as string, argv as list of string)` | Build a job from a name and an argv. |
| `pipeline.shell(name as string, cmdline as string)` | Build a job that runs a shell command line. |
| `pipeline.needing(j as Job, inputs as list of string)` | A copy of the job that declares the files it reads. |
| `pipeline.producing(j as Job, outputs as list of string)` | A copy of the job that declares the files it produces. |
| `pipeline.feeding(j as Job, text as string)` | A copy of the job with text to feed to its standard input. |
| `pipeline.tolerant(j as Job)` | A copy of the job whose failure is reported but not fatal. |
| `pipeline.needsRun(j as Job)` | Whether a job still needs to run. |
| `pipeline.run(j as Job)` | Run one job, skipping it when its outputs are already up to date. |
| `pipeline.runAll(jobs as list of Job)` | Run jobs one after another, stopping at the first real failure. |
| `pipeline.runParallel(jobs as list of Job, width as int)` | Run jobs concurrently, at most `width` at a time. |
| `pipeline.check(res as Result)` | Raise when a result represents a real failure. |
| `pipeline.checkAll(results as list of Result)` | Raise on the first real failure in a batch of results. |
| `pipeline.allOk(results as list of Result)` | Whether every result succeeded or was tolerated. |
| `pipeline.hasTool(name as string)` | Whether an executable is on `PATH`. |
| `pipeline.requireTools(names as list of string)` | Raise unless every named tool is available. |
| `pipeline.workDir(tag as string)` | A scratch directory that the caller is responsible for removing. |
| `pipeline.cleanUp(dir as string)` | Remove a directory and everything under it. |
| `pipeline.summary(results as list of Result)` | Render results as one status line per job. |


## `ngs`

The deck entry module: identity and host capabilities.


### Types

| type | kind | description |
| --- | --- | --- |
| `ngs.Capabilities` | struct | What this host can and cannot do with the deck. |


### Constants

| constant | type | value | description |
| --- | --- | --- | --- |
| `ngs.VERSION` | `string` | `"0.1.0"` | The deck's version. |
| `ngs.NAME` | `string` | `"@jennifer/ngs"` | The deck's canonical package name. |


### Functions

| function | description |
| --- | --- |
| `ngs.capabilities()` | What the current interpreter and host support. |
| `ngs.describe()` | A one-line description of the deck and what this host supports. |
| `ngs.requireFullHost()` | Raise unless the host can stream gzip and drive external tools. |


---

244 exported names across 12 modules.
