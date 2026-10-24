# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Shared primitives: error construction, Go-speed counting, tab-delimited
 * field access and sequence helpers.
 *
 * Every other module imports this one; it imports nothing from the deck.
 *
 * @module core
 * @since 0.1.0
 */
use strings;
use convert;
use os;
use fs;
use path;

/** Prefix of every error kind the deck raises: `ngs.<module>`. */
export def const ERROR_KIND as string init "ngs";

/** Field separator of every text format in the deck. */
export def const TAB as string init "\t";

/**
 * Zero run used to left-pad sort keys.
 *
 * Slicing the tail of `PAD_ZEROS + digits` costs one allocation, where
 * `strings.repeat` costs a length computation and a build; over the hundreds of
 * thousands of keys a sort of a real annotation needs, that is measurable.
 * 19 digits covers any 64-bit coordinate.
 */
def const PAD_ZEROS as string init "0000000000000000000";

/**
 * Raise a deck error.
 *
 * @param what {string} subsystem tag, e.g. "source" or "fastq"
 * @param message {string} human-readable detail
 * @throws {Error} kind "ngs.<what>", always
 */
export func fail(what as string, message as string) {
    throw Error{kind: ERROR_KIND + "." + $what, message: $message, file: "", line: 0, col: 0};
}

/**
 * Count non-overlapping occurrences of `needle` in `haystack`.
 *
 * Runs inside Go (`strings.replace` plus two `len`s): one native scan instead
 * of a per-character Jennifer loop. `qc` counts base composition and adapters
 * with it.
 *
 * @param haystack {string} text to scan
 * @param needle {string} substring to count
 * @return {int} number of occurrences, 0 when `needle` is empty
 */
export func countIn(haystack as string, needle as string) {
    if (len($needle) == 0) {
        return 0;
    }
    return (len($haystack) - len(strings.replace($haystack, $needle, ""))) // len($needle);
}

/**
 * Parse an integer field, reporting a malformed value as a deck error.
 *
 * `convert.toInt` raises a bare `runtime` error naming only the text; this
 * raises `ngs.<what>` naming the column. The record parsers convert columns
 * inline and wrap the whole record with `reraise` instead.
 *
 * @param what {string} subsystem tag, e.g. "bed"
 * @param label {string} what the column is, e.g. "start"
 * @param raw {string} the text to parse
 * @return {int} the parsed value
 * @throws {Error} kind "ngs.<what>" when the text is not an integer
 */
export func parseInt(what as string, label as string, raw as string) {
    try {
        return convert.toInt($raw);
    } catch (e) {
        fail($what, $label + " is not an integer: \"" + $raw + "\"");
    }
    return 0;
}

/**
 * Parse a float field, reporting a malformed value as a deck error.
 *
 * @param what {string} subsystem tag, e.g. "vcf"
 * @param label {string} what the column is, e.g. "qual"
 * @param raw {string} the text to parse
 * @return {float} the parsed value
 * @throws {Error} kind "ngs.<what>" when the text is not a number
 */
export func parseFloat(what as string, label as string, raw as string) {
    try {
        return convert.toFloat($raw);
    } catch (e) {
        fail($what, $label + " is not a number: \"" + $raw + "\"");
    }
    return 0.0;
}

/**
 * Re-raise a parse failure as a deck error naming the record it came from.
 *
 * Wrapping a whole record parse in one `try` and re-raising here costs a single
 * exception frame per record, where guarding each numeric column separately
 * cost one per column - measurably, on files with millions of records. An error
 * that already carries an `ngs.*` kind is re-thrown untouched, so the module
 * that diagnosed it keeps the credit.
 *
 * @param what {string} subsystem tag, e.g. "bed"
 * @param line {string} record being parsed, quoted in the message
 * @param e {Error} error that interrupted the parse
 * @throws {Error} always: `e` itself when its kind starts with "ngs.", else kind
 *   "ngs.<what>" with the record appended to the message
 */
export func reraise(what as string, line as string, e as Error) {
    if (strings.startsWith($e.kind, ERROR_KIND + ".")) {
        throw $e;
    }
    fail($what, $e.message + " in record: " + $line);
}

/**
 * Split a tab-delimited line into its fields.
 *
 * @param line {string} one record, newline already stripped
 * @return {list of string} the fields; never empty
 */
export func splitTab(line as string) {
    return strings.split($line, TAB);
}

/**
 * Field `i` of `fields`, or `fallback` when the column is absent.
 *
 * Optional trailing columns are normal in BED and SAM.
 *
 * @param fields {list of string} split record
 * @param i {int} zero-based column index
 * @param fallback {string} value for a missing column
 * @return {string} the column, or `fallback`
 */
export func field(fields as list of string, i as int, fallback as string) {
    if ($i < 0 or $i >= len($fields)) {
        return $fallback;
    }
    return $fields[$i];
}

/**
 * Integer field `i`, or `fallback` when absent, empty or a placeholder dot.
 *
 * @param what {string} subsystem tag for the error kind, e.g. "gff"
 * @param fields {list of string} split record
 * @param i {int} zero-based column index
 * @param fallback {int} value for a missing or "." column
 * @return {int} the parsed column, or `fallback`
 * @throws {Error} kind "ngs.<what>" when the column is present but not an integer
 */
export func intField(what as string, fields as list of string, i as int, fallback as int) {
    def raw as string init field($fields, $i, "");
    if (len($raw) == 0 or $raw == ".") {
        return $fallback;
    }
    return parseInt($what, "column " + convert.toString($i + 1), $raw);
}

/**
 * Float field `i`, or `fallback` when absent, empty or a placeholder dot.
 *
 * @param what {string} subsystem tag for the error kind, e.g. "vcf"
 * @param fields {list of string} split record
 * @param i {int} zero-based column index
 * @param fallback {float} value for a missing or "." column
 * @return {float} the parsed column, or `fallback`
 * @throws {Error} kind "ngs.<what>" when the column is present but not a number
 */
export func floatField(what as string, fields as list of string, i as int, fallback as float) {
    def raw as string init field($fields, $i, "");
    if (len($raw) == 0 or $raw == ".") {
        return $fallback;
    }
    return parseFloat($what, "column " + convert.toString($i + 1), $raw);
}

/**
 * Whether a record reader should skip a line.
 *
 * @param line {string} raw line
 * @param marker {string} comment prefix, e.g. "#" or "@"
 * @return {bool} true for a blank or whitespace-only line, or one starting with `marker`
 */
export func isSkippable(line as string, marker as string) {
    if (len(strings.trim($line)) == 0) {
        return true;
    }
    return strings.startsWith($line, $marker);
}

/**
 * Left-pad an integer with zeros to `width`, for building sortable string keys.
 *
 * Negative values are returned unpadded (genomic coordinates are non-negative),
 * so they sort before every padded value.
 *
 * @param n {int} value to pad
 * @param width {int} minimum number of digits
 * @return {string} the digits, zero-padded on the left
 */
export func padNum(n as int, width as int) {
    if ($n < 0) {
        return convert.toString($n);
    }
    def s as string init convert.toString($n);
    if (len($s) >= $width) {
        return $s;
    }
    if ($width > len(PAD_ZEROS)) {
        # Beyond the constant's reach, fall back to building the run. This is
        # not a path any coordinate takes - it exists so an exported function is
        # total rather than raising a substring error from its own internals.
        return strings.repeat("0", $width - len($s)) + $s;
    }
    def padded as string init PAD_ZEROS + $s;
    return strings.substring($padded, len($padded) - $width);
}

/**
 * Quote a string for safe inclusion in a POSIX `sh -c` command line.
 *
 * Wraps in single quotes and escapes embedded single quotes, so a path
 * containing spaces, `$`, `;` or quotes cannot break out into shell syntax.
 *
 * @param s {string} literal value to quote
 * @return {string} shell-safe token including its surrounding quotes
 */
export func shellQuote(s as string) {
    return "'" + strings.replace($s, "'", "'\\''") + "'";
}

/**
 * Reverse complement of a DNA sequence.
 *
 * Complements A, C, G, T, U (to A) and the IUPAC ambiguity codes in either
 * case; any other character, N included, is kept as is.
 *
 * @param seq {string} nucleotide sequence
 * @return {string} reverse complement, case preserved
 */
export func reverseComplement(seq as string) {
    def out as list of string init [];
    def cs as list of string init strings.chars($seq);
    def i as int init len($cs) - 1;
    while ($i >= 0) {
        $out[] = complementBase($cs[$i]);
        $i = $i - 1;
    }
    return strings.join($out, "");
}

/**
 * Complement of a single nucleotide character.
 *
 * @param b {string} one-character base
 * @return {string} complementary base, or `b` when it has none
 */
export func complementBase(b as string) {
    match ($b) {
        when "A" { return "T"; }
        when "C" { return "G"; }
        when "G" { return "C"; }
        when "T" { return "A"; }
        when "U" { return "A"; }
        when "a" { return "t"; }
        when "c" { return "g"; }
        when "g" { return "c"; }
        when "t" { return "a"; }
        when "u" { return "a"; }
        when "R" { return "Y"; }
        when "Y" { return "R"; }
        when "S" { return "S"; }
        when "W" { return "W"; }
        when "K" { return "M"; }
        when "M" { return "K"; }
        when "B" { return "V"; }
        when "V" { return "B"; }
        when "D" { return "H"; }
        when "H" { return "D"; }
        when "r" { return "y"; }
        when "y" { return "r"; }
        when "k" { return "m"; }
        when "m" { return "k"; }
        when "b" { return "v"; }
        when "v" { return "b"; }
        when "d" { return "h"; }
        when "h" { return "d"; }
        else { return $b; }
    }
}

/**
 * G+C count of a sequence, at Go speed.
 *
 * @param seq {string} nucleotide sequence
 * @return {int} number of G, C, g and c characters
 */
export func gcCount(seq as string) {
    return countIn($seq, "G") + countIn($seq, "C") + countIn($seq, "g") + countIn($seq, "c");
}

/**
 * N count of a sequence, at Go speed.
 *
 * @param seq {string} nucleotide sequence
 * @return {int} number of N and n characters
 */
export func nCount(seq as string) {
    return countIn($seq, "N") + countIn($seq, "n");
}

/**
 * Whether a program is on `PATH`.
 *
 * Walks `PATH` with `fs` rather than shelling out, so it works without the
 * `exec` capability. Checks only that a regular file of that name exists, not
 * that it is executable.
 *
 * @param name {string} program name, or a path containing a separator
 * @return {bool} true when a regular file of that name is found
 */
export func hasTool(name as string) {
    if (strings.contains($name, os.DIRSEP)) {
        return fs.isFile($name);
    }
    for (def dir in strings.split(os.getEnv("PATH"), os.PATHSEP)) {
        if (len($dir) == 0) {
            continue;
        }
        if (fs.isFile(path.join($dir, $name))) {
            return true;
        }
    }
    return false;
}

/**
 * Whether two paths resolve to the same file.
 *
 * Resolves symlinks, so `./reads.fq`, `reads.fq` and a link to it compare
 * equal. A path that does not exist is never the same file as anything, itself
 * included, so a not-yet-created output gives `false`.
 *
 * @param a {string} first path
 * @param b {string} second path
 * @return {bool} true when both exist and resolve to the same file
 */
export func sameFile(a as string, b as string) {
    if (not fs.exists($a) or not fs.exists($b)) {
        return false;
    }
    return fs.realpath($a) == fs.realpath($b);
}
