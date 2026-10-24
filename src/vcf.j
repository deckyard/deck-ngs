# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming VCF reading and writing.
 *
 * `open` captures the `##` meta lines and the `#CHROM` sample names on the
 * `Reader`; `next` yields only variant records.
 *
 * **Coordinates.** `Variant.pos` is VCF's 1-based POS; `toInterval` converts to
 * 0-based half-open, spanning the reference allele.
 *
 * @module vcf
 * @since 0.1.0
 */
use strings;
use convert;
use maps;
use lists;
import "./core.j" as core;
import "./source.j" as source;
import "./intervals.j" as intervals;

/** Placeholder for an absent column. */
export def const MISSING as string init ".";

/** Columns before the first sample column: the eight fixed ones plus FORMAT. */
export def const FIXED_COLUMNS as int init 9;

/**
 * One VCF variant record.
 *
 * @field chrom {string} reference sequence name
 * @field pos {int} 1-based position of the first reference base
 * @field id {string} variant identifier, "" when "."
 * @field ref {string} reference allele
 * @field alt {list of string} alternate alleles, empty when "."
 * @field qual {float} Phred-scaled quality, -1.0 when "."
 * @field filter {string} FILTER column verbatim, e.g. "PASS"
 * @field info {map of string to string} parsed INFO; a flag maps to ""
 * @field format {list of string} FORMAT keys, in column order
 * @field samples {list of string} per-sample columns verbatim
 */
export def struct Variant {
    chrom as string,
    pos as int,
    id as string,
    ref as string,
    alt as list of string,
    qual as float,
    filter as string,
    info as map of string to string,
    format as list of string,
    samples as list of string
};

/**
 * A streaming VCF cursor, carrying the file's header.
 *
 * @field src {source.Source} underlying line source
 * @field meta {list of string} `##` meta-information lines, in file order
 * @field samples {list of string} sample names from the `#CHROM` line
 */
export def struct Reader {
    src as source.Source,
    meta as list of string,
    samples as list of string
};

/**
 * The result of advancing a VCF cursor.
 *
 * @field Record {Variant} next variant, in `value`
 * @field Done {none} end of stream
 */
export def enum Next {
    Record{value as Variant},
    Done
};

/**
 * Open a VCF file for streaming, capturing its header.
 *
 * The header comes from a separate short pass over the front of the file
 * (`headerOf`); the record stream starts at byte zero and `next` skips the `#`
 * lines. A gzipped file is decompressed from the start twice.
 *
 * @param path {string} `.vcf` or gzipped equivalent (plain gzip or BGZF)
 * @return {Reader} an open cursor with `meta` and `samples` populated
 * @throws {Error} kind "ngs.source" when the file is missing
 */
export func open(path as string) {
    def head as list of string init headerOf($path);
    return Reader{src: source.open($path), meta: metaOf($head), samples: samplesOf($head)};
}

/**
 * Open standard input as a VCF stream.
 *
 * `next` skips the header, but `meta` and `samples` stay empty: reading them
 * eagerly would consume the first variant, and a stream cannot be re-read. Use
 * `open` on a path when you need sample names. See `source.openStdin` for when
 * opening fails.
 *
 * @return {Reader} an open cursor over stdin, with an empty header
 */
export func openStdin() {
    def none as list of string init [];
    return Reader{src: source.openStdin(), meta: $none, samples: $none};
}

/**
 * Advance the cursor to the next variant record.
 *
 * Skips blank lines and `#` lines.
 *
 * @param r {Reader} an open cursor
 * @return {Next} `Record` with the variant, or `Done` at end of stream
 * @throws {Error} as `parse`
 */
export func next(r as Reader) {
    while (source.hasNext($r.src)) {
        def line as string init source.nextLine($r.src);
        if (not core.isSkippable($line, "#")) {
            return Next.Record{value: parse($line)};
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
 * Read the header lines of a VCF file and close it. Reads only the front of
 * the file; the decompressor's `SIGPIPE` exit is ignored.
 *
 * @param path {string} `.vcf` or gzipped equivalent
 * @return {list of string} every leading `#` line, in file order
 * @throws {Error} kind "ngs.source" when the file is missing
 */
export func headerOf(path as string) {
    def s as source.Source init source.open($path);
    errdefer quietClose($s);
    def lines as list of string init [];
    while (source.hasNext($s)) {
        def line as string init source.nextLine($s);
        if (not strings.startsWith($line, "#")) {
            break;
        }
        $lines[] = $line;
    }
    source.close($s);
    return $lines;
}

/**
 * Parse one VCF variant line.
 *
 * @param line {string} one record, newline already stripped
 * @return {Variant} the parsed variant
 * @throws {Error} kind "ngs.vcf" on fewer than eight columns or a non-numeric
 *   POS or QUAL
 */
export func parse(line as string) {
    try {
        return parseRow($line);
    } catch (e) {
        return core.reraise("vcf", $line, $e);
    }
}

/**
 * Parse a VCF record; the caller wraps it for error reporting.
 *
 * @param line {string} one record, newline already stripped
 * @return {Variant} the parsed variant
 * @throws {Error} kind "ngs.vcf" on fewer than eight columns; unwrapped
 *   conversion errors otherwise
 */
func parseRow(line as string) {
    def f as list of string init core.splitTab($line);
    if (len($f) < 8) {
        def got as string init convert.toString(len($f));
        core.fail("vcf", "VCF record needs 8+ columns, got " + $got + ": " + $line);
    }
    def id as string init $f[2];
    if ($id == MISSING) {
        $id = "";
    }
    def alts as list of string init [];
    if ($f[4] != MISSING) {
        if (strings.indexOf($f[4], ",") < 0) {
            $alts[] = $f[4];
        } else {
            $alts = strings.split($f[4], ",");
        }
    }
    def keys as list of string init [];
    if (len($f) > 8 and $f[8] != MISSING) {
        $keys = strings.split($f[8], ":");
    }
    def cols as list of string init [];
    if (len($f) > FIXED_COLUMNS) {
        $cols = lists.slice($f, FIXED_COLUMNS);
    }
    def quality as float init -1.0;
    if ($f[5] != MISSING and len($f[5]) > 0) {
        $quality = convert.toFloat($f[5]);
    }
    return Variant{
        chrom: $f[0],
        pos: convert.toInt($f[1]),
        id: $id,
        ref: $f[3],
        alt: $alts,
        qual: $quality,
        filter: $f[6],
        info: parseInfo($f[7]),
        format: $keys,
        samples: $cols
    };
}

/**
 * Parse a VCF INFO column into a map.
 *
 * A `key=value` entry maps to the text after the first `=`; a flag maps to "",
 * so test a flag with `hasInfo`. Values are not split on `,`.
 *
 * @param text {string} raw INFO column
 * @return {map of string to string} keys to values in column order, empty for "." or ""
 */
export func parseInfo(text as string) {
    def out as map of string to string init {};
    if ($text == MISSING or len($text) == 0) {
        return $out;
    }
    for (def entry in strings.split($text, ";")) {
        if (len($entry) == 0) {
            continue;
        }
        def cut as int init strings.indexOf($entry, "=");
        if ($cut < 0) {
            $out[$entry] = "";
            continue;
        }
        $out[strings.substring($entry, 0, $cut)] = strings.substring($entry, $cut + 1);
    }
    return $out;
}

/**
 * Read an INFO field, or a fallback when it is absent.
 *
 * @param v {Variant} the variant
 * @param key {string} INFO key
 * @param fallback {string} value for a missing key
 * @return {string} the value ("" for a flag), or `fallback`
 */
export func info(v as Variant, key as string, fallback as string) {
    if (not maps.has($v.info, $key)) {
        return $fallback;
    }
    return $v.info[$key];
}

/**
 * Whether an INFO key or flag is present.
 *
 * @param v {Variant} the variant
 * @param key {string} INFO key
 * @return {bool} true when present
 */
export func hasInfo(v as Variant, key as string) {
    return maps.has($v.info, $key);
}

/**
 * Whether the variant passed every filter.
 *
 * `PASS`, `.` and an empty FILTER all count as passing.
 *
 * @param v {Variant} the variant
 * @return {bool} true when passing
 */
export func passed(v as Variant) {
    return $v.filter == "PASS" or $v.filter == MISSING or len($v.filter) == 0;
}

/**
 * Whether the variant is a single-nucleotide substitution.
 *
 * True only when the reference and every alternate allele are one of the
 * upper-case bases A, C, G, T or N. A padding or symbolic reference such as
 * `-` is not an SNV, even when the alternate is a plain base; lower-case
 * alleles do not count.
 *
 * @param v {Variant} the variant
 * @return {bool} true for an SNV
 */
export func isSnv(v as Variant) {
    if (len($v.ref) != 1 or len($v.alt) == 0 or not isBase($v.ref)) {
        return false;
    }
    for (def a in $v.alt) {
        if (len($a) != 1 or not isBase($a)) {
            return false;
        }
    }
    return true;
}

/**
 * Whether the variant is an insertion or deletion.
 *
 * True when any alternate allele is a plain sequence that differs in length
 * from the reference. Symbolic (`<DEL>`) and breakend (`A[chr2:123[`) alleles
 * do not count (see `isSequence`); the `*` allele is not excluded.
 *
 * @param v {Variant} the variant
 * @return {bool} true for an indel
 */
export func isIndel(v as Variant) {
    for (def a in $v.alt) {
        if (len($a) != len($v.ref) and isSequence($a)) {
            return true;
        }
    }
    return false;
}

/**
 * The per-sample value of one FORMAT key.
 *
 * @param v {Variant} the variant
 * @param sampleIndex {int} zero-based index into the reader's `samples`
 * @param key {string} FORMAT key, e.g. "GT" or "DP"
 * @param fallback {string} value when the sample, key or trailing field is absent
 * @return {string} the value, or `fallback`
 */
export func sampleValue(v as Variant, sampleIndex as int, key as string, fallback as string) {
    if ($sampleIndex < 0 or $sampleIndex >= len($v.samples)) {
        return $fallback;
    }
    def slot as int init indexOfKey($v.format, $key);
    if ($slot < 0) {
        return $fallback;
    }
    def parts as list of string init strings.split($v.samples[$sampleIndex], ":");
    if ($slot >= len($parts)) {
        return $fallback;
    }
    return $parts[$slot];
}

/**
 * The genotype call of one sample.
 *
 * @param v {Variant} the variant
 * @param sampleIndex {int} zero-based index into the reader's `samples`
 * @return {string} the GT value as written, e.g. "0/1" or "1|0"; "./." when absent
 */
export func genotype(v as Variant, sampleIndex as int) {
    return sampleValue($v, $sampleIndex, "GT", "./.");
}

/**
 * The reference span of a variant as a 0-based half-open interval:
 * `pos - 1` to that plus the length of `ref`.
 *
 * @param v {Variant} the variant
 * @return {Interval} the range, named by the id, QUAL as score (0.0 when
 *   missing), strand "."
 * @throws {Error} kind "ngs.intervals" when `pos` is below 1
 */
export func toInterval(v as Variant) {
    def start as int init $v.pos - 1;
    def score as float init $v.qual;
    if ($score < 0.0) {
        $score = 0.0;
    }
    return intervals.named($v.chrom, $start, $start + len($v.ref), $v.id, $score, ".");
}

/**
 * Render a variant as a VCF line. FORMAT is written only when `format` is not
 * empty; QUAL goes through `convert.toString`, so 60 is written "60.0".
 *
 * @param v {Variant} the variant
 * @return {string} one tab-delimited record, no trailing newline
 */
export func format(v as Variant) {
    def id as string init $v.id;
    if (len($id) == 0) {
        $id = MISSING;
    }
    def alt as string init MISSING;
    if (len($v.alt) > 0) {
        $alt = strings.join($v.alt, ",");
    }
    def qual as string init MISSING;
    if ($v.qual >= 0.0) {
        $qual = convert.toString($v.qual);
    }
    def columns as list of string init [
        $v.chrom,
        core.padNum($v.pos, 1),
        $id,
        $v.ref,
        $alt,
        $qual,
        $v.filter,
        formatInfo($v.info)
    ];
    if (len($v.format) > 0) {
        $columns[] = strings.join($v.format, ":");
    }
    for (def col in $v.samples) {
        $columns[] = $col;
    }
    return strings.join($columns, core.TAB);
}

/**
 * Render an INFO map as column text; an entry with an empty value is written
 * as a flag.
 *
 * @param entries {map of string to string} INFO keys to values, in output order
 * @return {string} the INFO column, "." when empty
 */
export func formatInfo(entries as map of string to string) {
    def parts as list of string init [];
    for (def key in $entries) {
        if (len($entries[$key]) == 0) {
            $parts[] = $key;
        } else {
            $parts[] = $key + "=" + $entries[$key];
        }
    }
    if (len($parts) == 0) {
        return MISSING;
    }
    return strings.join($parts, ";");
}

/**
 * The `##` meta lines of a captured header.
 *
 * @param head {list of string} header lines from `headerOf`
 * @return {list of string} the `##` lines
 */
func metaOf(head as list of string) {
    def out as list of string init [];
    for (def line in $head) {
        if (strings.startsWith($line, "##")) {
            $out[] = $line;
        }
    }
    return $out;
}

/**
 * The sample names declared on the `#CHROM` line of a captured header.
 *
 * @param head {list of string} header lines from `headerOf`
 * @return {list of string} sample names, empty when the file has no genotypes
 */
func samplesOf(head as list of string) {
    def out as list of string init [];
    for (def line in $head) {
        if (not strings.startsWith($line, "#CHROM")) {
            continue;
        }
        def f as list of string init core.splitTab($line);
        if (len($f) <= FIXED_COLUMNS) {
            continue;
        }
        for (def i in FIXED_COLUMNS..len($f)) {
            $out[] = $f[$i];
        }
    }
    return $out;
}

/**
 * Position of a key in the FORMAT list.
 *
 * @param keys {list of string} FORMAT keys
 * @param key {string} key to find
 * @return {int} zero-based index, -1 when absent
 */
func indexOfKey(keys as list of string, key as string) {
    for (def i in 0..len($keys)) {
        if ($keys[$i] == $key) {
            return $i;
        }
    }
    return -1;
}

/**
 * Whether a one-character allele is an upper-case nucleotide code.
 *
 * @param a {string} the allele
 * @return {bool} true for A, C, G, T and N
 */
func isBase(a as string) {
    return $a == "A" or $a == "C" or $a == "G" or $a == "T" or $a == "N";
}

/**
 * Whether an allele is a plain sequence rather than a symbolic or breakend
 * allele such as `<DEL>` or `A[chr2:123[`. Only `<`, `[` and `]` are checked.
 *
 * @param a {string} the allele
 * @return {bool} false also for an empty allele
 */
func isSequence(a as string) {
    if (len($a) == 0) {
        return false;
    }
    return not strings.startsWith($a, "<") and strings.indexOf($a, "[") < 0 and
        strings.indexOf($a, "]") < 0;
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
