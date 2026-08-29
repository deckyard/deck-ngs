# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming VCF reading and writing.
 *
 * The `##` meta lines and the `#CHROM` column header are captured when the
 * reader is opened, so sample names are available up front and `next` yields
 * only variant records.
 *
 * **Coordinates.** VCF `POS` is 1-based and inclusive; `toInterval` converts to
 * the 0-based half-open convention used by `intervals`, spanning the reference
 * allele.
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

/** Placeholder used for an absent column. */
export def const MISSING as string init ".";

/** Number of fixed columns before the per-sample genotype columns. */
export def const FIXED_COLUMNS as int init 9;

/**
 * One VCF variant record.
 *
 * @field chrom {string} reference sequence name
 * @field pos {int} 1-based position of the first reference base
 * @field id {string} variant identifier, "" when absent
 * @field ref {string} the reference allele
 * @field alt {list of string} alternate alleles, empty when "."
 * @field qual {float} phred-scaled variant quality, -1.0 when absent
 * @field filter {string} the FILTER column verbatim, e.g. "PASS"
 * @field info {map of string to string} parsed INFO; a flag maps to ""
 * @field format {list of string} the FORMAT keys, in column order
 * @field samples {list of string} the per-sample columns verbatim
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
 * @field src {source.Source} the underlying line source
 * @field meta {list of string} the `##` meta-information lines, in file order
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
 * @field Record {Variant} the next parsed variant, in `value`
 * @field Done {none} the stream is exhausted
 */
export def enum Next {
    Record{value as Variant},
    Done
};

/**
 * Open a VCF file for streaming, capturing its header.
 *
 * The header is read by a short separate pass over the front of the file, so
 * sample names are known before the first record while the record stream still
 * starts clean; `next` skips the `#` lines as it walks.
 *
 * @param path {string} `.vcf` or gzipped equivalent
 * @return {Reader} an open cursor with `meta` and `samples` populated
 * @throws {Error} when the file is missing
 */
export func open(path as string) {
    def head as list of string init headerOf($path);
    return Reader{src: source.open($path), meta: metaOf($head), samples: samplesOf($head)};
}

/**
 * Open standard input as a VCF stream.
 *
 * `next` skips the header, but `meta` and `samples` are left empty: a stream
 * cannot be re-read, and reading the header eagerly would consume the first
 * variant with no way to put it back. Use `open` on a path when you need
 * sample names.
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
 * @param r {Reader} an open cursor
 * @return {Next} `Record` with the variant, or `Done` at end of stream
 * @throws {Error} on a record with fewer than eight columns
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
 * Read just the header lines of a VCF file.
 *
 * @param path {string} `.vcf` or gzipped equivalent
 * @return {list of string} every leading `#` line, in file order
 * @throws {Error} when the file is missing
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
 * @throws {Error} when fewer than eight columns are present
 */
export func parse(line as string) {
    try {
        return parseRow($line);
    } catch (e) {
        return core.reraise("vcf", $line, $e);
    }
}

/**
 * Parse a VCF record, assuming the caller has wrapped it for error reporting.
 *
 * @param line {string} one record, newline already stripped
 * @return {Variant} the parsed variant
 * @throws {Error} when fewer than eight columns are present
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
 * A `key=value` entry maps to its value; a bare flag maps to the empty string,
 * so use `maps.has` to test for a flag's presence.
 *
 * @param text {string} the raw INFO column
 * @return {map of string to string}
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
 * @param fallback {string} value to return when the key is missing
 * @return {string}
 */
export func info(v as Variant, key as string, fallback as string) {
    if (not maps.has($v.info, $key)) {
        return $fallback;
    }
    return $v.info[$key];
}

/**
 * Whether an INFO flag is present.
 *
 * @param v {Variant} the variant
 * @param key {string} INFO key
 * @return {bool}
 */
export func hasInfo(v as Variant, key as string) {
    return maps.has($v.info, $key);
}

/**
 * Whether the variant passed every filter.
 *
 * Treats both `PASS` and `.` as passing, matching common caller output.
 *
 * @param v {Variant} the variant
 * @return {bool}
 */
export func passed(v as Variant) {
    return $v.filter == "PASS" or $v.filter == MISSING or len($v.filter) == 0;
}

/**
 * Whether the variant is a single-nucleotide substitution.
 *
 * True only when the reference and every alternate allele are a single
 * nucleotide. A padding or symbolic reference such as `-` is not an SNV, even
 * when the alternate is a plain base.
 *
 * @param v {Variant} the variant
 * @return {bool}
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
 * True when any allele differs in length from the reference.
 *
 * @param v {Variant} the variant
 * @return {bool}
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
 * @param fallback {string} value to return when the key or column is absent
 * @return {string}
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
 * @return {string} the GT value, e.g. "0/1", or "./." when unavailable
 */
export func genotype(v as Variant, sampleIndex as int) {
    return sampleValue($v, $sampleIndex, "GT", "./.");
}

/**
 * The reference span of a variant, as a half-open interval.
 *
 * Converts VCF's 1-based POS to the 0-based convention; the interval covers the
 * reference allele.
 *
 * @param v {Variant} the variant
 * @return {Interval} the reference range, named by the variant id
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
 * Render a variant back to a VCF record line.
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
 * Render an INFO map back to its column text.
 *
 * @param entries {map of string to string} INFO keys to values
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
 * @return {list of string} only the `##` lines
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
 * @param keys {list of string} the FORMAT keys
 * @param key {string} the key to find
 * @return {int} zero-based index, or -1 when absent
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
 * Whether a one-character allele is a nucleotide code.
 *
 * @param a {string} the allele
 * @return {bool}
 */
func isBase(a as string) {
    return $a == "A" or $a == "C" or $a == "G" or $a == "T" or $a == "N";
}

/**
 * Whether an allele is a plain nucleotide sequence rather than a symbolic or
 * breakend allele such as `<DEL>` or `A[chr2:123[`.
 *
 * @param a {string} the allele
 * @return {bool}
 */
func isSequence(a as string) {
    if (len($a) == 0) {
        return false;
    }
    return not strings.startsWith($a, "<") and strings.indexOf($a, "[") < 0 and
        strings.indexOf($a, "]") < 0;
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
