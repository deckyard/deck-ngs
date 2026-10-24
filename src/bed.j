# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming BED reading and writing. BED and `intervals.Interval` are both
 * 0-based half-open, so a record parses straight into an interval.
 *
 * `track`, `browser` and `#` lines may sit between records, and skipping one
 * can consume the last line of the file, which a `hasNext` peek cannot
 * express. The cursor returns a `Next` enum instead:
 *
 * ```jennifer
 * def r as bed.Reader init bed.open("peaks.bed");
 * def going as bool init true;
 * while ($going) {
 *     def step as bed.Next init bed.next($r);
 *     match ($step) {
 *         when Feature(f) { io.printf("%s\n", intervals.toBed3($f.value)); }
 *         when Done { $going = false; }
 *     }
 * }
 * bed.closeStrict($r);
 * ```
 *
 * `readAll` loads a whole file, for feature sets that fit in memory.
 *
 * @module bed
 * @since 0.1.0
 */
use strings;
use convert;
import "./core.j" as core;
import "./source.j" as source;
import "./intervals.j" as intervals;

/** Placeholder read and written for an absent name, score or strand. */
export def const MISSING as string init ".";

/**
 * A streaming BED cursor.
 *
 * @field src {source.Source} underlying line source
 */
export def struct Reader {
    src as source.Source
};

/**
 * The result of advancing a BED cursor.
 *
 * @field Feature {Interval} next record, in `value`
 * @field Done {none} end of stream
 */
export def enum Next {
    Feature{value as intervals.Interval},
    Done
};

/**
 * Open a BED file for streaming, decompressing gzip transparently.
 *
 * @param path {string} `.bed` or gzipped equivalent
 * @return {Reader} an open cursor
 * @throws {Error} kind "ngs.source" when the file is missing
 */
export func open(path as string) {
    return Reader{src: source.open($path)};
}

/**
 * Open standard input as a BED stream. See `source.openStdin` for when that
 * fails.
 *
 * @return {Reader} an open cursor over stdin
 */
export func openStdin() {
    return Reader{src: source.openStdin()};
}

/**
 * Advance the cursor to the next feature, skipping blank, comment, `track` and
 * `browser` lines.
 *
 * @param r {Reader} an open cursor
 * @return {Next} `Feature` with the record, or `Done` at end of stream
 * @throws {Error} as `parse`
 */
export func next(r as Reader) {
    while (source.hasNext($r.src)) {
        def line as string init source.nextLine($r.src);
        if (not isDecoration($line)) {
            return Next.Feature{value: parse($line)};
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
 * Parse one BED line into an interval.
 *
 * Reads the first six columns and ignores the rest (thickStart, blocks, ...).
 * A missing or "." name becomes "", a missing or "." score 0.0, a missing
 * strand ".".
 *
 * @param line {string} one BED record, newline already stripped
 * @return {Interval} the parsed feature
 * @throws {Error} kind "ngs.bed" on fewer than three columns or a non-numeric
 *   coordinate or score; kind "ngs.intervals" on negative or inverted
 *   coordinates
 */
export func parse(line as string) {
    try {
        return parseRow($line);
    } catch (e) {
        return core.reraise("bed", $line, $e);
    }
}

/**
 * Parse a BED record; the caller wraps it for error reporting.
 *
 * Reads columns inline rather than through `core.field` and `core.floatField`:
 * those are two to three extra calls per record, which is the difference
 * between a fast and a slow pass over a million features.
 *
 * @param line {string} one BED record, newline already stripped
 * @return {Interval} the parsed feature
 * @throws {Error} kind "ngs.bed" on fewer than three columns; unwrapped
 *   conversion errors otherwise
 */
func parseRow(line as string) {
    def f as list of string init core.splitTab($line);
    if (len($f) < 3) {
        def got as string init convert.toString(len($f));
        core.fail("bed", "BED record needs at least 3 columns, got " + $got + ": " + $line);
    }
    def name as string init "";
    if (len($f) > 3 and $f[3] != MISSING) {
        $name = $f[3];
    }
    def score as float init 0.0;
    if (len($f) > 4 and len($f[4]) > 0 and $f[4] != MISSING) {
        $score = convert.toFloat($f[4]);
    }
    def strand as string init MISSING;
    if (len($f) > 5) {
        $strand = $f[5];
    }
    return intervals.named(
        $f[0],
        convert.toInt($f[1]),
        convert.toInt($f[2]),
        $name,
        $score,
        $strand);
}

/**
 * Render an interval as a BED6 line. An empty name or strand is written ".";
 * a whole-number score is written without a decimal part.
 *
 * @param iv {Interval} the feature
 * @return {string} tab-delimited BED6, no trailing newline
 */
export func format(iv as intervals.Interval) {
    def name as string init $iv.name;
    if (len($name) == 0) {
        $name = MISSING;
    }
    def strand as string init $iv.strand;
    if (len($strand) == 0) {
        $strand = MISSING;
    }
    return strings.join(
        [
            $iv.chrom,
            core.padNum($iv.start, 1),
            core.padNum($iv.end, 1),
            $name,
            formatScore($iv.score),
            $strand
        ],
        core.TAB);
}

/**
 * Read a whole BED file into memory.
 *
 * Memory grows with the file; for files that do not fit, use `next`.
 *
 * @param path {string} `.bed` or gzipped equivalent
 * @return {list of Interval} every feature, in file order
 * @throws {Error} kind "ngs.source" when the file is missing or the gzip stream
 *   is truncated; as `parse` on a malformed record
 */
export func readAll(path as string) {
    def r as Reader init open($path);
    errdefer quietClose($r);
    def out as list of intervals.Interval init [];
    # Reads lines directly rather than driving the `Next` cursor: the enum wraps
    # each feature in a payload that is copied into the binder and copied again
    # out of it, which doubles the struct copies on the one path whose whole job
    # is to build a list of them.
    while (source.hasNext($r.src)) {
        def line as string init source.nextLine($r.src);
        if (not isDecoration($line)) {
            $out[] = parse($line);
        }
    }
    closeStrict($r);
    return $out;
}

/**
 * Write intervals to a BED file, gzipping when the name ends in `.gz`.
 *
 * @param path {string} file to create or overwrite
 * @param ivs {list of Interval} features to write, in order
 * @throws {Error} kind "ngs.bed" when the compressor exits non-zero; as
 *   `source.create` when the target cannot be opened
 */
export func writeAll(path as string, ivs as list of intervals.Interval) {
    def k as source.Sink init source.create($path);
    errdefer quietCloseSink($k);
    for (def iv in $ivs) {
        source.writeLine($k, format($iv));
    }
    def code as int init source.closeSink($k);
    if ($code != 0) {
        core.fail(
            "bed",
            "compressor for " + $path + " exited with code " + convert.toString($code));
    }
}

/**
 * Write one interval to an open sink as a `format` line.
 *
 * @param k {source.Sink} an open sink
 * @param iv {Interval} feature to write
 */
export func write(k as source.Sink, iv as intervals.Interval) {
    source.writeLine($k, format($iv));
}

/**
 * Whether a line is a BED decoration rather than a record.
 *
 * @param line {string} raw line
 * @return {bool} true for blank, `#`, `track` and `browser` lines
 */
func isDecoration(line as string) {
    if (len($line) == 0) {
        return true;
    }
    def first as string init strings.substring($line, 0, 1);
    if ($first == "#") {
        return true;
    }
    if ($first == "t") {
        return strings.startsWith($line, "track");
    }
    if ($first == "b") {
        return strings.startsWith($line, "browser");
    }
    if ($first == " " or $first == core.TAB) {
        return len(strings.trim($line)) == 0;
    }
    return false;
}

/**
 * Render a score, dropping the decimal part when it is a whole number.
 *
 * @param score {float} feature score
 * @return {string} score text
 */
func formatScore(score as float) {
    def rounded as int init convert.toInt($score);
    if ($rounded == $score) {
        return convert.toString($rounded);
    }
    return convert.toString($score);
}

/**
 * Close a cursor on a path that is already failing, without raising.
 *
 * @param r {Reader} cursor to close
 * @return {int} the close result, -1 when closing itself failed
 */
func quietClose(r as Reader) {
    try {
        return close($r);
    } catch (e) {
        return -1;
    }
}

/**
 * Close a sink on a path that is already failing, without raising.
 *
 * @param k {source.Sink} sink to close
 * @return {int} the close result, -1 when closing itself failed
 */
func quietCloseSink(k as source.Sink) {
    try {
        return source.closeSink($k);
    } catch (e) {
        return -1;
    }
}
