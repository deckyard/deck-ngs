# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming SAM reading and writing.
 *
 * SAM only. BAM and CRAM need BGZF (blocked gzip with random access), which
 * `compress` does not provide; pipe them through `samtools view` and read the
 * SAM output here.
 *
 * `open` captures the `@` header block on the `Reader`; `next` yields only
 * alignment records.
 *
 * **Coordinates.** `Record.pos` is 1-based as in the file; `toInterval`
 * converts to 0-based half-open and takes the reference end from the CIGAR.
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

/** FLAG 0x1: template has multiple segments (paired). */
export def const FLAG_PAIRED as int init 0x1;

/** FLAG 0x2: each segment properly aligned, according to the aligner. */
export def const FLAG_PROPER_PAIR as int init 0x2;

/** FLAG 0x4: segment unmapped. */
export def const FLAG_UNMAPPED as int init 0x4;

/** FLAG 0x8: next segment (mate) unmapped. */
export def const FLAG_MATE_UNMAPPED as int init 0x8;

/** FLAG 0x10: segment reverse complemented. */
export def const FLAG_REVERSE as int init 0x10;

/** FLAG 0x20: mate reverse complemented. */
export def const FLAG_MATE_REVERSE as int init 0x20;

/** FLAG 0x40: first segment in the template. */
export def const FLAG_FIRST as int init 0x40;

/** FLAG 0x80: last segment in the template. */
export def const FLAG_LAST as int init 0x80;

/** FLAG 0x100: secondary alignment. */
export def const FLAG_SECONDARY as int init 0x100;

/** FLAG 0x200: failed platform or vendor quality checks. */
export def const FLAG_QC_FAIL as int init 0x200;

/** FLAG 0x400: PCR or optical duplicate. */
export def const FLAG_DUPLICATE as int init 0x400;

/** FLAG 0x800: supplementary alignment. */
export def const FLAG_SUPPLEMENTARY as int init 0x800;

/** Placeholder for an absent string column. */
export def const MISSING as string init "*";

/**
 * One SAM alignment record.
 *
 * @field qname {string} query template name
 * @field flag {int} bitwise FLAG
 * @field rname {string} reference sequence name, "*" when unmapped
 * @field pos {int} 1-based leftmost mapping position, 0 when unmapped
 * @field mapq {int} mapping quality, 255 when unavailable
 * @field cigar {string} raw CIGAR, "*" when unavailable
 * @field rnext {string} reference name of the mate, "=" when the same
 * @field pnext {int} 1-based position of the mate
 * @field tlen {int} observed template length, signed
 * @field seq {string} segment sequence, "*" when absent
 * @field qual {string} raw ASCII quality string, "*" when absent
 * @field tags {list of string} optional `TAG:TYPE:VALUE` columns, verbatim
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
 * @field length {int} bases the operation covers
 * @field op {string} one-character code (M, I, D, N, S, H, P, =, X)
 */
export def struct CigarOp {
    length as int,
    op as string
};

/**
 * A streaming SAM cursor, carrying the file's header block.
 *
 * @field src {source.Source} underlying line source
 * @field header {list of string} `@` header lines, in file order
 */
export def struct Reader {
    src as source.Source,
    header as list of string
};

/**
 * The result of advancing a SAM cursor.
 *
 *   `Alignment`  the next record, a `Record`, in `value`
 *   `Done`       the end of the stream
 */
export def enum Next {
    Alignment{value as Record},
    Done
};

/**
 * Open a SAM file for streaming, capturing its header block.
 *
 * The header comes from a separate short pass over the front of the file
 * (`headerOf`); the record stream starts at byte zero and `next` skips the `@`
 * lines. A gzipped file is decompressed from the start twice.
 *
 * @param path {string} `.sam` or gzipped equivalent
 * @return {Reader} an open cursor with `header` populated
 * @throws {Error} kind "ngs.source" when the file is missing
 */
export func open(path as string) {
    def head as list of string init headerOf($path);
    return Reader{src: source.open($path), header: $head};
}

/**
 * Open standard input as a SAM stream.
 *
 * For `samtools view -h` output. `next` skips the header block, but `header`
 * stays empty: reading it eagerly would consume the first alignment, and a
 * stream cannot be re-read. Use `open` on a path when you need
 * `referenceSizes`. See `source.openStdin` for when opening fails.
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
 * Skips blank lines and `@` lines.
 *
 * @param r {Reader} an open cursor
 * @return {Next} `Alignment` with the record, or `Done` at end of stream
 * @throws {Error} as `parse`
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
 * Close the cursor and release its decompressor. See `source.close`.
 *
 * @param r {Reader} cursor to close
 * @return {int} the decompressor's exit code, 0 when uncompressed
 */
export func close(r as Reader) {
    return source.close($r.src);
}

/**
 * Close a fully-consumed cursor, raising if its decompressor failed.
 *
 * @param r {Reader} cursor to close, already read to exhaustion
 * @throws {Error} kind "ngs.source" on a truncated or corrupt gzip stream
 */
export func closeStrict(r as Reader) {
    source.closeStrict($r.src);
}

/**
 * Parse one SAM alignment line.
 *
 * @param line {string} one record, newline already stripped
 * @return {Record} the parsed alignment
 * @throws {Error} kind "ngs.sam" on fewer than eleven columns or a non-numeric
 *   FLAG, POS, MAPQ, PNEXT or TLEN
 */
export func parse(line as string) {
    try {
        return parseRow($line);
    } catch (e) {
        return core.reraise("sam", $line, $e);
    }
}

/**
 * Parse a SAM record; the caller wraps it for error reporting.
 *
 * Converts numeric columns directly rather than through a guarded helper per
 * column: one exception frame for the record replaces five, which is the bulk
 * of the per-record cost on a large alignment file.
 *
 * @param line {string} one record, newline already stripped
 * @return {Record} the parsed alignment
 * @throws {Error} kind "ngs.sam" on fewer than eleven columns; unwrapped
 *   conversion errors otherwise
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
 * Render an alignment as a SAM line, tags appended verbatim.
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
 * @param bit {int} a `FLAG_*` constant, or several OR-ed together
 * @return {bool} true when any of the bits is set
 */
export func hasFlag(rec as Record, bit as int) {
    return ($rec.flag & $bit) != 0;
}

/**
 * Whether the record is unmapped (0x4).
 *
 * @param rec {Record} the alignment
 * @return {bool} true when unmapped
 */
export func isUnmapped(rec as Record) {
    return hasFlag($rec, FLAG_UNMAPPED);
}

/**
 * Whether the record is reverse complemented (0x10).
 *
 * @param rec {Record} the alignment
 * @return {bool} true on the reverse strand
 */
export func isReverse(rec as Record) {
    return hasFlag($rec, FLAG_REVERSE);
}

/**
 * Whether the record is a secondary alignment (0x100).
 *
 * @param rec {Record} the alignment
 * @return {bool} true when secondary
 */
export func isSecondary(rec as Record) {
    return hasFlag($rec, FLAG_SECONDARY);
}

/**
 * Whether the record is a supplementary alignment (0x800).
 *
 * @param rec {Record} the alignment
 * @return {bool} true when supplementary
 */
export func isSupplementary(rec as Record) {
    return hasFlag($rec, FLAG_SUPPLEMENTARY);
}

/**
 * Whether the record is a PCR or optical duplicate (0x400).
 *
 * @param rec {Record} the alignment
 * @return {bool} true when flagged as a duplicate
 */
export func isDuplicate(rec as Record) {
    return hasFlag($rec, FLAG_DUPLICATE);
}

/**
 * Whether the record is a primary alignment: none of 0x4, 0x100, 0x800 set.
 * The usual filter for counting reads once; duplicates and QC failures still
 * pass.
 *
 * @param rec {Record} the alignment
 * @return {bool} true when primary
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
 * Operation codes are not validated; any non-digit character is taken as one.
 *
 * @param cigar {string} raw CIGAR, "*" for none
 * @return {list of CigarOp} operations in order, empty for "*" or ""
 * @throws {Error} kind "ngs.sam" on an operation without a length or a trailing
 *   length without an operation
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
 * Counts M, D, N, = and X; I, S, H and P do not advance the reference.
 *
 * @param cigar {string} raw CIGAR
 * @return {int} reference span in bases, 0 for "*"
 * @throws {Error} as `parseCigar`
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
 * @param cigar {string} raw CIGAR
 * @return {int} query span in bases, 0 for "*"
 * @throws {Error} as `parseCigar`
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
 * The reference range an alignment covers, as a 0-based half-open interval:
 * `pos - 1` to that plus `referenceLength(cigar)`. When the CIGAR spans no
 * reference bases ("*"), the length of `seq` is used instead.
 *
 * @param rec {Record} a mapped alignment
 * @return {Interval} the range, named after the read, MAPQ as score, strand
 *   from `strandOf`
 * @throws {Error} kind "ngs.sam" when the record is unmapped or `rname` is "*",
 *   or as `parseCigar`
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
 * @param tagName {string} two-character tag, e.g. "NM"
 * @param fallback {string} value for an absent tag
 * @return {string} the value as text, `TAG:TYPE:` removed, or `fallback`
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
 * Parses the `@SQ SN:name LN:length` lines into the size map
 * `intervals.complement` takes. Empty for a reader from `openStdin`.
 *
 * @param r {Reader} a cursor from `open`
 * @return {map of string to int} reference name to length, in header order
 * @throws {Error} kind "ngs.sam" when an LN value is not an integer
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
 * Read only the `@` header block of a SAM file.
 *
 * Reads up to the first line not starting with `@` and closes, so it touches
 * only the front of the file. Uses `close`, not `closeStrict`: abandoning a
 * gzip stream early makes the decompressor exit on `SIGPIPE`, which is
 * expected here.
 *
 * @param path {string} `.sam` or gzipped equivalent
 * @return {list of string} the header lines, in file order
 * @throws {Error} kind "ngs.source" when the file is missing
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
 * @param line {string} header line
 * @param prefix {string} field prefix including its colon
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
 * @param op {string} operation code
 * @return {bool} true for M, D, N, = and X
 */
func consumesReference(op as string) {
    return $op == "M" or $op == "D" or $op == "N" or $op == "=" or $op == "X";
}

/**
 * Whether a CIGAR operation consumes query bases.
 *
 * @param op {string} operation code
 * @return {bool} true for M, I, S, = and X
 */
func consumesQuery(op as string) {
    return $op == "M" or $op == "I" or $op == "S" or $op == "=" or $op == "X";
}

/**
 * Whether a one-character string is an ASCII digit.
 *
 * @param ch {string} one character
 * @return {bool} true for 0-9
 */
func isDigit(ch as string) {
    return $ch >= "0" and $ch <= "9";
}

/**
 * Close a source on a path that is already failing, without raising.
 *
 * @param s {source.Source} source to close
 * @return {int} the close result, -1 when closing itself failed
 */
func quietClose(s as source.Source) {
    try {
        return source.close($s);
    } catch (e) {
        return -1;
    }
}
