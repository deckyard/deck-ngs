# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

/**
 * Streaming FASTQ reading, writing and per-record accessors. A `Reader` is a
 * pull cursor over a `source.Source`, so gzipped input streams at constant
 * memory:
 *
 * ```jennifer
 * def r as fastq.Reader init fastq.open("sample.fastq.gz");
 * while (fastq.hasNext($r)) {
 *     def rec as fastq.Read init fastq.next($r);
 *     io.printf("%s %d\n", $rec.id, fastq.length($rec));
 * }
 * fastq.closeStrict($r);
 * ```
 *
 * Records are four-line FASTQ; multi-line (wrapped) records are not read.
 * Quality strings stay raw ASCII. Decoding them to Phred scores costs one
 * Jennifer step per base, so `qualities`, `meanQuality` and `minQuality` are
 * the expensive accessors; for whole-file statistics use `qc`.
 *
 * @module fastq
 * @since 0.1.0
 */
use strings;
use convert;
use fs;
use channel;
use regex;
import "./core.j" as core;
import "./source.j" as source;

/** Sanger / Illumina 1.8+ quality offset, the default everywhere in the deck. */
export def const PHRED33 as int init 33;

/** Legacy Illumina 1.3-1.7 quality offset. */
export def const PHRED64 as int init 64;

/**
 * One FASTQ record.
 *
 * @field id {string} read name: the header up to the first space, without `@`
 * @field desc {string} header text after the first space, "" when none
 * @field seq {string} nucleotide sequence
 * @field qual {string} raw ASCII quality string, same length as `seq`, every
 *   byte in `!`..`~` (33..126)
 */
export def struct Read {
    id as string,
    desc as string,
    seq as string,
    qual as string
};

/**
 * A streaming FASTQ cursor. Build it with `open` or `openStdin`.
 *
 * @field src {source.Source} underlying line source
 * @field ahead {channel of string} the header line `hasNext` read ahead, for
 *   `next` to take; a handle, so copies of the reader share it
 */
export def struct Reader {
    src as source.Source,
    ahead as channel of string
};

# Any byte outside printable ASCII `!`..`~`, the FASTQ quality alphabet.
def const NOT_QUALITY as string init "[^!-~]";

/**
 * Open a FASTQ file for streaming, decompressing gzip transparently.
 *
 * @param path {string} `.fastq`, `.fq` or gzipped equivalent
 * @return {Reader} an open cursor; close it with `close` or `closeStrict`
 * @throws {Error} kind "ngs.source" when the file is missing, or is gzipped and
 *   the host lacks `exec`
 */
export func open(path as string) {
    return Reader{src: source.open($path), ahead: channel.make(1)};
}

/**
 * Open standard input as a FASTQ stream. See `source.openStdin` for when that
 * fails.
 *
 * @return {Reader} an open cursor over stdin
 */
export func openStdin() {
    return Reader{src: source.openStdin(), ahead: channel.make(1)};
}

/**
 * Whether another record is available.
 *
 * Reads ahead to the next non-blank line, which `next` then takes as the
 * header, so blank lines at the end of the file (or between records) are
 * skipped and do not count as a record.
 *
 * @param r {Reader} an open cursor
 * @return {bool} false at end of stream, or when only blank lines remain
 */
export func hasNext(r as Reader) {
    if (channel.len($r.ahead) > 0) {
        return true;
    }
    while (not fs.eof($r.src.file)) {
        def line as string init fs.readLine($r.src.file);
        if (strings.startsWith($line, "@") or len(strings.trim($line)) > 0) {
            channel.send($r.ahead, $line);
            return true;
        }
    }
    return false;
}

/**
 * Read the next record.
 *
 * Blank lines before the header are skipped.
 *
 * @param r {Reader} an open cursor
 * @return {Read} the next record
 * @throws {Error} kind "ngs.fastq" at end of stream, on a truncated record, a
 *   malformed header or separator, a sequence/quality length mismatch or a
 *   quality byte outside `!`..`~`, naming the header
 */
export func next(r as Reader) {
    # `fs` directly rather than `source.hasNext`/`nextLine`: in this loop a
    # call into another module costs more than the read itself.
    def header as string init "";
    try {
        if (channel.len($r.ahead) > 0) {
            $header = channel.recv($r.ahead);
        } else {
            $header = fs.readLine($r.src.file);
            while (not strings.startsWith($header, "@") and len(strings.trim($header)) == 0) {
                $header = fs.readLine($r.src.file);
            }
        }
        # Reading the three remaining lines under one exception frame replaces a
        # `hasNext` probe between each of them. A truncated record then surfaces
        # as the read error it is, re-raised with the record's header for
        # context - and the probes were a quarter of the cost of the hottest
        # loop in the deck.
        def seq as string init fs.readLine($r.src.file);
        def plus as string init fs.readLine($r.src.file);
        def qual as string init fs.readLine($r.src.file);
        return parse($header, $seq, $plus, $qual);
    } catch (e) {
        return core.reraise("fastq", $header, $e);
    }
}

/**
 * Close the reader and release its decompressor. See `source.close`.
 *
 * @param r {Reader} cursor to close
 * @return {int} the decompressor's exit code, 0 when uncompressed
 */
export func close(r as Reader) {
    return source.close($r.src);
}

/**
 * Close a fully-consumed reader, raising if its decompressor failed.
 *
 * @param r {Reader} cursor to close, already read to exhaustion
 * @throws {Error} kind "ngs.source" on a truncated or corrupt gzip stream
 */
export func closeStrict(r as Reader) {
    source.closeStrict($r.src);
}

/**
 * Build a record from its four raw lines, validating the FASTQ shape.
 *
 * @param header {string} `@name description` line
 * @param seq {string} sequence line
 * @param plus {string} separator line, starting with `+`
 * @param qual {string} quality line
 * @return {Read} the parsed record
 * @throws {Error} kind "ngs.fastq" when a marker is wrong, a quality byte is
 *   outside `!`..`~` (33..126), or sequence and quality differ in length in
 *   bytes
 */
export func parse(header as string, seq as string, plus as string, qual as string) {
    if (not strings.startsWith($header, "@")) {
        core.fail("fastq", "record header must start with '@', got: " + $header);
    }
    if (not strings.startsWith($plus, "+")) {
        core.fail("fastq", "record separator must start with '+' in read: " + $header);
    }
    # Inline rather than `checkQuality`: one call fewer per record.
    if (regex.matches(NOT_QUALITY, $qual)) {
        core.fail("fastq", "quality outside '!'..'~' (ASCII 33..126) in read: " + $header);
    }
    # The quality is ASCII now, so its rune count is its byte count.
    if (len($seq) != len($qual) or len(convert.bytesFromString($seq, "utf-8")) != len($qual)) {
        core.fail("fastq", "sequence and quality differ in length in read: " + $header);
    }
    def rest as string init strings.substring($header, 1);
    def cut as int init strings.indexOf($rest, " ");
    if ($cut < 0) {
        return Read{id: $rest, desc: "", seq: $seq, qual: $qual};
    }
    return Read{
        id: strings.substring($rest, 0, $cut),
        desc: strings.substring($rest, $cut + 1),
        seq: $seq,
        qual: $qual
    };
}

/**
 * Render a record as four-line FASTQ. The separator is always a bare `+`.
 *
 * @param rec {Read} the record
 * @return {string} four lines, newline-terminated
 * @throws {Error} kind "ngs.fastq" when the record would not read back as
 *   itself: whitespace in `id`, CR or LF in `desc`, tab, CR or LF in `seq`, a
 *   quality byte outside `!`..`~`, or `seq` and `qual` of different byte length
 */
export func format(rec as Read) {
    checkHeader($rec);
    checkSequence($rec);
    checkQuality($rec.qual, $rec.id);
    if (len($rec.seq) != len($rec.qual) or
        len(convert.bytesFromString($rec.seq, "utf-8")) != len($rec.qual)) {
        core.fail("fastq", "sequence and quality differ in length in read: " + $rec.id);
    }
    return "@" + headerOf($rec) + "\n" + $rec.seq + "\n+\n" + $rec.qual + "\n";
}

/**
 * Render a record as two-line FASTA, quality discarded.
 *
 * @param rec {Read} the record
 * @return {string} two lines, newline-terminated
 * @throws {Error} kind "ngs.fastq" on whitespace in `id`, CR or LF in `desc`, or
 *   tab, CR or LF in `seq`
 */
export func toFasta(rec as Read) {
    checkHeader($rec);
    checkSequence($rec);
    return ">" + headerOf($rec) + "\n" + $rec.seq + "\n";
}

/**
 * Reject an id or description that would break the header line on write.
 * A tab in `desc` is allowed: `samtools fastq -T` writes SAM tags there.
 *
 * @param rec {Read} the record
 * @throws {Error} kind "ngs.fastq" on whitespace in `id` or CR/LF in `desc`
 */
func checkHeader(rec as Read) {
    if (regex.matches("\\s", $rec.id)) {
        core.fail("fastq", "read id contains whitespace: " + $rec.id);
    }
    if (strings.indexOf($rec.desc, "\n") >= 0 or strings.indexOf($rec.desc, "\r") >= 0) {
        core.fail("fastq", "read description contains a line break in read: " + $rec.id);
    }
}

/**
 * Reject a sequence that would break the record's line structure.
 *
 * @param rec {Read} the record
 * @throws {Error} kind "ngs.fastq" on a tab, CR or LF in `seq`
 */
func checkSequence(rec as Read) {
    if (regex.matches("[\\t\\r\\n]", $rec.seq)) {
        core.fail("fastq", "sequence contains a tab or line break in read: " + $rec.id);
    }
}

/**
 * Reject a quality string with a byte outside the FASTQ alphabet `!`..`~`.
 * One native regex scan; a non-ASCII rune fails it too, so the string's rune
 * count is its byte count afterwards.
 *
 * @param qual {string} raw quality string
 * @param label {string} header or id for the message
 * @throws {Error} kind "ngs.fastq" on any byte outside 33..126
 */
func checkQuality(qual as string, label as string) {
    if (regex.matches(NOT_QUALITY, $qual)) {
        core.fail("fastq", "quality outside '!'..'~' (ASCII 33..126) in read: " + $label);
    }
}

/**
 * The header text of a record: id, plus a space and the description when
 * there is one.
 *
 * @param rec {Read} the record
 * @return {string} header without its leading marker
 */
export func headerOf(rec as Read) {
    if (len($rec.desc) == 0) {
        return $rec.id;
    }
    return $rec.id + " " + $rec.desc;
}

/**
 * Write a record to a sink as `format` renders it.
 *
 * @param k {source.Sink} an open sink
 * @param rec {Read} record to write
 */
export func write(k as source.Sink, rec as Read) {
    source.writeText($k, format($rec));
}

/**
 * Read length.
 *
 * @param rec {Read} the record
 * @return {int} length in bases
 */
export func length(rec as Read) {
    return len($rec.seq);
}

/**
 * Decode the quality string to Phred scores. Costs one Jennifer step per base.
 *
 * @param rec {Read} the record
 * @param offset {int} `PHRED33` or `PHRED64`
 * @return {list of int} one Phred score per base
 */
export func qualities(rec as Read, offset as int) {
    def raw as bytes init convert.bytesFromString($rec.qual, "utf-8");
    def out as list of int init [];
    for (def i in 0..len($raw)) {
        $out[] = $raw[$i] - $offset;
    }
    return $out;
}

/**
 * Mean Phred quality of a record.
 *
 * The arithmetic mean of the Phred scores, as FastQC and fastp report per
 * read, not the mean of error probabilities.
 *
 * @param rec {Read} the record
 * @param offset {int} `PHRED33` or `PHRED64`
 * @return {float} mean quality, 0.0 for an empty read
 */
export func meanQuality(rec as Read, offset as int) {
    if (len($rec.qual) == 0) {
        return 0.0;
    }
    def raw as bytes init convert.bytesFromString($rec.qual, "utf-8");
    def total as int init 0;
    for (def i in 0..len($raw)) {
        $total = $total + ($raw[$i] - $offset);
    }
    return $total / len($raw);
}

/**
 * Lowest Phred score in a record.
 *
 * @param rec {Read} the record
 * @param offset {int} `PHRED33` or `PHRED64`
 * @return {int} minimum quality, -1 for an empty read
 */
export func minQuality(rec as Read, offset as int) {
    if (len($rec.qual) == 0) {
        return -1;
    }
    def raw as bytes init convert.bytesFromString($rec.qual, "utf-8");
    def lowest as int init $raw[0] - $offset;
    for (def i in 0..len($raw)) {
        def q as int init $raw[$i] - $offset;
        if ($q < $lowest) {
            $lowest = $q;
        }
    }
    return $lowest;
}

/**
 * G+C content of a record, counted at Go speed.
 *
 * @param rec {Read} the record
 * @return {float} percentage in [0, 100], 0.0 for an empty read
 */
export func gcPercent(rec as Read) {
    if (len($rec.seq) == 0) {
        return 0.0;
    }
    return 100.0 * core.gcCount($rec.seq) / len($rec.seq);
}

/**
 * Number of N bases in a record, counted at Go speed.
 *
 * @param rec {Read} the record
 * @return {int} count of N and n
 */
export func ambiguousCount(rec as Read) {
    return core.nCount($rec.seq);
}

/**
 * A sub-range of a record, slicing sequence and quality together.
 *
 * @param rec {Read} the record
 * @param start {int} zero-based inclusive start
 * @param end {int} zero-based exclusive end
 * @return {Read} a record over the half-open range, id and desc preserved
 * @throws {Error} kind "ngs.fastq" when the range is outside the read or inverted
 */
export func subseq(rec as Read, start as int, end as int) {
    if ($start < 0 or $end > len($rec.seq) or $start > $end) {
        core.fail("fastq", "subseq range out of bounds for read: " + $rec.id);
    }
    return Read{
        id: $rec.id,
        desc: $rec.desc,
        seq: strings.substring($rec.seq, $start, $end),
        qual: strings.substring($rec.qual, $start, $end)
    };
}

/**
 * Reverse complement of a record, quality string reversed to match.
 *
 * @param rec {Read} the record
 * @return {Read} the reverse-complemented record
 */
export func reverseComplement(rec as Read) {
    def qs as list of string init strings.chars($rec.qual);
    def flipped as list of string init [];
    def i as int init len($qs) - 1;
    while ($i >= 0) {
        $flipped[] = $qs[$i];
        $i = $i - 1;
    }
    return Read{
        id: $rec.id,
        desc: $rec.desc,
        seq: core.reverseComplement($rec.seq),
        qual: strings.join($flipped, "")
    };
}
