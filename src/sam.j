# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming SAM reading and writing.
 *
 * SAM only - BAM and CRAM are deliberately out of scope for this deck: they need
 * BGZF (blocked gzip with random access), which `compress` does not provide, and
 * their value lies in the alignment compute this deck does not attempt. Pipe
 * through `samtools view` and read its SAM output here.
 *
 * The `@`-prefixed header block is consumed when the reader is opened and kept
 * on the `Reader`, so `next` only ever yields alignment records.
 *
 * **Coordinates.** SAM `pos` is 1-based; `toInterval` converts to the 0-based
 * half-open convention and uses the CIGAR string to compute the reference end.
 *
 * @module sam
 * @since 0.1.0
 */
use strings;
use convert;
use lists;
import "./core.j" as core;
import "./source.j" as source;
import "./intervals.j" as intervals;

/** Read is paired in sequencing. */
export def const FLAG_PAIRED as int init 0x1;

/** Each segment is properly aligned according to the aligner. */
export def const FLAG_PROPER_PAIR as int init 0x2;

/** This segment is unmapped. */
export def const FLAG_UNMAPPED as int init 0x4;

/** The mate segment is unmapped. */
export def const FLAG_MATE_UNMAPPED as int init 0x8;

/** This segment is on the reverse strand. */
export def const FLAG_REVERSE as int init 0x10;

/** The mate segment is on the reverse strand. */
export def const FLAG_MATE_REVERSE as int init 0x20;

/** This is the first segment of the template. */
export def const FLAG_FIRST as int init 0x40;

/** This is the last segment of the template. */
export def const FLAG_LAST as int init 0x80;

/** This is a secondary alignment. */
export def const FLAG_SECONDARY as int init 0x100;

/** The read failed platform or vendor quality checks. */
export def const FLAG_QC_FAIL as int init 0x200;

/** The read is a PCR or optical duplicate. */
export def const FLAG_DUPLICATE as int init 0x400;

/** This is a supplementary alignment. */
export def const FLAG_SUPPLEMENTARY as int init 0x800;

/** Placeholder used for an absent column. */
export def const MISSING as string init "*";

/**
 * One SAM alignment record.
 *
 * @field qname {string} query template name
 * @field flag {int} bitwise FLAG field
 * @field rname {string} reference sequence name, "*" when unmapped
 * @field pos {int} 1-based leftmost mapping position, 0 when unmapped
 * @field mapq {int} mapping quality
 * @field cigar {string} the raw CIGAR string, "*" when unavailable
 * @field rnext {string} reference name of the mate
 * @field pnext {int} 1-based position of the mate
 * @field tlen {int} observed template length, signed
 * @field seq {string} the segment sequence, "*" when absent
 * @field qual {string} raw ASCII quality string, "*" when absent
 * @field tags {list of string} the optional `TAG:TYPE:VALUE` columns, verbatim
 */
export def struct Record {
    qname as string,
    flag as int,
    rname as string,
    pos as int,
    mapq as int,
    cigar as string,
    rnext as string,
    pnext as int,
    tlen as int,
    seq as string,
    qual as string,
    tags as list of string
};

/**
 * One CIGAR operation.
 *
 * @field length {int} number of bases the operation covers
 * @field op {string} the single-character operation code (M, I, D, N, S, H, P, =, X)
 */
export def struct CigarOp {
    length as int,
    op as string
};

/**
 * A streaming SAM cursor, carrying the file's header block.
 *
 * @field src {source.Source} the underlying line source
 * @field header {list of string} the `@`-prefixed header lines, in file order
 */
export def struct Reader {
    src as source.Source,
    header as list of string
};

/**
 * The result of advancing a SAM cursor.
 *
 * @field Alignment {Record} the next parsed record, in `value`
 * @field Done {none} the stream is exhausted
 */
export def enum Next {
    Alignment{value as Record},
    Done
};

/**
 * Open a SAM file for streaming, capturing its header block.
 *
 * The header is read by a short separate pass over the front of the file, so the
 * record stream itself starts clean at byte zero; `next` then skips the `@`
 * lines as it walks. Headers are small, so the second pass costs little.
 *
 * @param path {string} `.sam` or gzipped equivalent
 * @return {Reader} an open cursor with `header` populated
 * @throws {Error} when the file is missing
 */
export func open(path as string) {
    def head as list of string init headerOf($path);
    return Reader{src: source.open($path), header: $head};
}

/**
 * Open standard input as a SAM stream.
 *
 * The natural way to read `samtools view -h` output. `next` skips the header
 * block, but `header` is left empty: a stream cannot be re-read, and reading the
 * header eagerly would consume the first alignment with no way to put it back.
 * Use `open` on a path when you need `referenceSizes`.
 *
 * @return {Reader} an open cursor over stdin, with an empty `header`
 */
export func openStdin() {
    def head as list of string init [];
    return Reader{src: source.openStdin(), header: $head};
}

/**
 * Advance the cursor to the next alignment record.
 *
 * @param r {Reader} an open cursor
 * @return {Next} `Alignment` with the record, or `Done` at end of stream
 * @throws {Error} on a record with fewer than eleven columns
 */
export func next(r as Reader) {
    while (source.hasNext($r.src)) {
        def line as string init source.nextLine($r.src);
        if (not core.isSkippable($line, "@")) {
            return Next.Alignment{value: parse($line)};
        }
    }
    return Next.Done;
}

/**
 * Close the cursor and release its decompressor.
 *
 * @param r {Reader} the cursor to close
 * @return {int} the decompressor's exit code, 0 when uncompressed
 */
export func close(r as Reader) {
    return source.close($r.src);
}

/**
 * Close a fully-consumed cursor, raising if its decompressor failed.
 *
 * @param r {Reader} the cursor to close
 * @throws {Error} on a truncated or corrupt gzip stream
 */
export func closeStrict(r as Reader) {
    source.closeStrict($r.src);
}

/**
 * Parse one SAM alignment line.
 *
 * @param line {string} one record, newline already stripped
 * @return {Record} the parsed alignment
 * @throws {Error} when fewer than eleven columns are present
 */
export func parse(line as string) {
    try {
        return parseRow($line);
    } catch (e) {
        return core.reraise("sam", $line, $e);
    }
}

/**
 * Parse a SAM record, assuming the caller has wrapped it for error reporting.
 *
 * Converts numeric columns directly rather than through a guarded helper per
 * column: one exception frame for the record replaces five, which is the bulk
 * of the per-record cost on a large alignment file.
 *
 * @param line {string} one record, newline already stripped
 * @return {Record} the parsed alignment
 * @throws {Error} when fewer than eleven columns are present
 */
func parseRow(line as string) {
    def f as list of string init core.splitTab($line);
    if (len($f) < 11) {
        def got as string init convert.toString(len($f));
        core.fail("sam", "SAM record needs 11+ columns, got " + $got + ": " + $line);
    }
    def tags as list of string init [];
    if (len($f) > 11) {
        $tags = lists.slice($f, 11);
    }
    return Record{
        qname: $f[0],
        flag: convert.toInt($f[1]),
        rname: $f[2],
        pos: convert.toInt($f[3]),
        mapq: convert.toInt($f[4]),
        cigar: $f[5],
        rnext: $f[6],
        pnext: convert.toInt($f[7]),
        tlen: convert.toInt($f[8]),
        seq: $f[9],
        qual: $f[10],
        tags: $tags
    };
}

/**
 * Render an alignment back to a SAM record line.
 *
 * @param rec {Record} the alignment
 * @return {string} one tab-delimited record, no trailing newline
 */
export func format(rec as Record) {
    def columns as list of string init [
        $rec.qname,
        core.padNum($rec.flag, 1),
        $rec.rname,
        core.padNum($rec.pos, 1),
        core.padNum($rec.mapq, 1),
        $rec.cigar,
        $rec.rnext,
        core.padNum($rec.pnext, 1),
        convert.toString($rec.tlen),
        $rec.seq,
        $rec.qual
    ];
    for (def tag in $rec.tags) {
        $columns[] = $tag;
    }
    return strings.join($columns, core.TAB);
}

/**
 * Whether a FLAG bit is set on a record.
 *
 * @param rec {Record} the alignment
 * @param bit {int} one of the `FLAG_*` constants
 * @return {bool}
 */
export func hasFlag(rec as Record, bit as int) {
    return ($rec.flag & $bit) != 0;
}

/**
 * Whether the record is unmapped.
 *
 * @param rec {Record} the alignment
 * @return {bool}
 */
export func isUnmapped(rec as Record) {
    return hasFlag($rec, FLAG_UNMAPPED);
}

/**
 * Whether the record is mapped to the reverse strand.
 *
 * @param rec {Record} the alignment
 * @return {bool}
 */
export func isReverse(rec as Record) {
    return hasFlag($rec, FLAG_REVERSE);
}

/**
 * Whether the record is a secondary alignment.
 *
 * @param rec {Record} the alignment
 * @return {bool}
 */
export func isSecondary(rec as Record) {
    return hasFlag($rec, FLAG_SECONDARY);
}

/**
 * Whether the record is a supplementary alignment.
 *
 * @param rec {Record} the alignment
 * @return {bool}
 */
export func isSupplementary(rec as Record) {
    return hasFlag($rec, FLAG_SUPPLEMENTARY);
}

/**
 * Whether the record is flagged as a PCR or optical duplicate.
 *
 * @param rec {Record} the alignment
 * @return {bool}
 */
export func isDuplicate(rec as Record) {
    return hasFlag($rec, FLAG_DUPLICATE);
}

/**
 * Whether the record is a primary alignment: mapped, not secondary, not
 * supplementary.
 *
 * The usual filter for counting reads once.
 *
 * @param rec {Record} the alignment
 * @return {bool}
 */
export func isPrimary(rec as Record) {
    # One mask instead of three calls into `hasFlag`; this runs once per record
    # in every alignment filter anyone writes.
    return ($rec.flag & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY)) == 0;
}

/**
 * The strand a record aligns to.
 *
 * @param rec {Record} the alignment
 * @return {string} "-" for reverse, "+" for forward, "." when unmapped
 */
export func strandOf(rec as Record) {
    if (($rec.flag & FLAG_UNMAPPED) != 0) {
        return ".";
    }
    if (($rec.flag & FLAG_REVERSE) != 0) {
        return "-";
    }
    return "+";
}

/**
 * Split a CIGAR string into its operations.
 *
 * @param cigar {string} the raw CIGAR, "*" for none
 * @return {list of CigarOp} operations in order, empty for "*"
 * @throws {Error} on a malformed CIGAR string
 */
export func parseCigar(cigar as string) {
    def out as list of CigarOp init [];
    if ($cigar == MISSING or len($cigar) == 0) {
        return $out;
    }
    def digits as string init "";
    for (def ch in strings.chars($cigar)) {
        if (isDigit($ch)) {
            $digits = $digits + $ch;
            continue;
        }
        if (len($digits) == 0) {
            core.fail("sam", "CIGAR operation without a length: " + $cigar);
        }
        $out[] = CigarOp{length: core.parseInt("sam", "CIGAR length", $digits), op: $ch};
        $digits = "";
    }
    if (len($digits) > 0) {
        core.fail("sam", "CIGAR ends with a length and no operation: " + $cigar);
    }
    return $out;
}

/**
 * Number of reference bases a CIGAR string consumes.
 *
 * Counts M, D, N, = and X; insertions, soft clips, hard clips and padding do
 * not advance the reference.
 *
 * @param cigar {string} the raw CIGAR
 * @return {int} reference span in bases
 */
export func referenceLength(cigar as string) {
    def total as int init 0;
    for (def step in parseCigar($cigar)) {
        if (consumesReference($step.op)) {
            $total = $total + $step.length;
        }
    }
    return $total;
}

/**
 * Number of query bases a CIGAR string consumes.
 *
 * Counts M, I, S, = and X.
 *
 * @param cigar {string} the raw CIGAR
 * @return {int} query span in bases
 */
export func queryLength(cigar as string) {
    def total as int init 0;
    for (def step in parseCigar($cigar)) {
        if (consumesQuery($step.op)) {
            $total = $total + $step.length;
        }
    }
    return $total;
}

/**
 * The reference interval an alignment covers, as a half-open interval.
 *
 * Converts SAM's 1-based `pos` to the 0-based convention and derives the end
 * from the CIGAR string.
 *
 * @param rec {Record} a mapped alignment
 * @return {Interval} the covered reference range, named after the read
 * @throws {Error} when the record is unmapped
 */
export func toInterval(rec as Record) {
    if (isUnmapped($rec) or $rec.rname == MISSING) {
        core.fail("sam", "cannot place an unmapped read on the reference: " + $rec.qname);
    }
    def span as int init referenceLength($rec.cigar);
    if ($span == 0) {
        $span = len($rec.seq);
    }
    def start as int init $rec.pos - 1;
    def score as float init convert.toFloat($rec.mapq);
    return intervals.named($rec.rname, $start, $start + $span, $rec.qname, $score, strandOf($rec));
}

/**
 * Read the value of an optional tag.
 *
 * @param rec {Record} the alignment
 * @param tagName {string} the two-character tag name, e.g. "NM"
 * @param fallback {string} value to return when the tag is absent
 * @return {string} the tag's value with its `TAG:TYPE:` prefix removed
 */
export func tagValue(rec as Record, tagName as string, fallback as string) {
    def prefix as string init $tagName + ":";
    for (def entry in $rec.tags) {
        if (strings.startsWith($entry, $prefix)) {
            def rest as string init strings.substring($entry, len($prefix));
            def cut as int init strings.indexOf($rest, ":");
            if ($cut < 0) {
                return $rest;
            }
            return strings.substring($rest, $cut + 1);
        }
    }
    return $fallback;
}

/**
 * The reference sequences declared in the header, as chromosome sizes.
 *
 * Parses the `@SQ SN:name LN:length` lines, giving the chromosome-size map that
 * `intervals.complement` needs.
 *
 * @param r {Reader} an open cursor
 * @return {map of string to int} reference name to length
 */
export func referenceSizes(r as Reader) {
    def out as map of string to int init {};
    for (def line in $r.header) {
        if (not strings.startsWith($line, "@SQ")) {
            continue;
        }
        def sn as string init headerTag($line, "SN:");
        def ln as string init headerTag($line, "LN:");
        if (len($sn) > 0 and len($ln) > 0) {
            $out[$sn] = core.parseInt("sam", "@SQ LN", $ln);
        }
    }
    return $out;
}

/**
 * Read just the `@` header block of a SAM file.
 *
 * Opens the file, reads until the first alignment line and closes again, so it
 * touches only the front of the file. Uses `close` rather than `closeStrict`
 * because abandoning a gzip stream early makes the decompressor exit on
 * `SIGPIPE`, which is expected here.
 *
 * @param path {string} `.sam` or gzipped equivalent
 * @return {list of string} the header lines, in file order
 * @throws {Error} when the file is missing
 */
export func headerOf(path as string) {
    def s as source.Source init source.open($path);
    errdefer quietClose($s);
    def lines as list of string init [];
    while (source.hasNext($s)) {
        def line as string init source.nextLine($s);
        if (not strings.startsWith($line, "@")) {
            break;
        }
        $lines[] = $line;
    }
    source.close($s);
    return $lines;
}

/**
 * Read the value of a header field such as `SN:` from an `@SQ` line.
 *
 * @param line {string} the header line
 * @param prefix {string} the field prefix including its colon
 * @return {string} the field value, "" when absent
 */
func headerTag(line as string, prefix as string) {
    for (def part in core.splitTab($line)) {
        if (strings.startsWith($part, $prefix)) {
            return strings.substring($part, len($prefix));
        }
    }
    return "";
}

/**
 * Whether a CIGAR operation advances the reference position.
 *
 * @param op {string} the operation code
 * @return {bool}
 */
func consumesReference(op as string) {
    return $op == "M" or $op == "D" or $op == "N" or $op == "=" or $op == "X";
}

/**
 * Whether a CIGAR operation consumes query bases.
 *
 * @param op {string} the operation code
 * @return {bool}
 */
func consumesQuery(op as string) {
    return $op == "M" or $op == "I" or $op == "S" or $op == "=" or $op == "X";
}

/**
 * Whether a one-character string is an ASCII digit.
 *
 * @param ch {string} one character
 * @return {bool}
 */
func isDigit(ch as string) {
    return $ch >= "0" and $ch <= "9";
}

/**
 * Close a source on a path that is already failing.
 *
 * Cleanup that runs while an error is propagating must not raise one of its own.
 *
 * @param s {source.Source} the source to close
 * @return {int} the close result, or -1 when closing itself failed
 */
func quietClose(s as source.Source) {
    try {
        return source.close($s);
    } catch (e) {
        return -1;
    }
}
