# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming GFF3 and GTF reading and writing.
 *
 * GFF and GTF share nine tab-delimited columns and differ only in how the
 * ninth encodes attributes: GFF3 writes `key=value;key=value` with
 * percent-escaping, GTF writes `key "value"; key "value";`. `parse` detects the
 * style per record, so a mixed directory needs no flag.
 *
 * **Coordinates.** GFF and GTF are 1-based and inclusive; BED and the interval
 * engine are 0-based and half-open. `toInterval` converts, and is the only
 * correct way to move a feature into `intervals` - subtracting one by hand at
 * call sites is how off-by-one bugs get into annotation pipelines.
 *
 * @module gff
 * @since 0.1.0
 */
use strings;
use convert;
use maps;
use encoding;
import "./core.j" as core;
import "./source.j" as source;
import "./intervals.j" as intervals;

/** Placeholder used for an absent column. */
export def const MISSING as string init ".";

/** Attribute style: `key=value;` pairs, percent-escaped. */
export def const STYLE_GFF3 as string init "gff3";

/** Attribute style: `key "value";` pairs. */
export def const STYLE_GTF as string init "gtf";

/**
 * One GFF3 or GTF feature.
 *
 * @field seqid {string} reference sequence name
 * @field origin {string} the `source` column: what produced the feature
 * @field kind {string} the `type` column: gene, exon, CDS, ...
 * @field start {int} 1-based inclusive start, as written in the file
 * @field end {int} 1-based inclusive end, as written in the file
 * @field score {float} feature score, 0.0 when absent
 * @field strand {string} "+", "-" or "."
 * @field phase {int} CDS reading frame 0-2, or -1 when absent
 * @field attributes {map of string to string} the parsed ninth column
 */
export def struct Feature {
    seqid as string,
    origin as string,
    kind as string,
    start as int,
    end as int,
    score as float,
    strand as string,
    phase as int,
    attributes as map of string to string
};

/**
 * A streaming GFF/GTF cursor.
 *
 * @field src {source.Source} the underlying line source
 */
export def struct Reader {
    src as source.Source
};

/**
 * The result of advancing a GFF cursor.
 *
 * @field Record {Feature} the next parsed feature, in `value`
 * @field Done {none} the stream is exhausted
 */
export def enum Next {
    Record{value as Feature},
    Done
};

/**
 * Open a GFF3 or GTF file for streaming, decompressing gzip transparently.
 *
 * @param path {string} `.gff`, `.gff3`, `.gtf` or gzipped equivalent
 * @return {Reader} an open cursor
 * @throws {Error} when the file is missing
 */
export func open(path as string) {
    return Reader{src: source.open($path)};
}

/**
 * Open standard input as a GFF/GTF stream.
 *
 * @return {Reader} an open cursor over stdin
 */
export func openStdin() {
    return Reader{src: source.openStdin()};
}

/**
 * Advance the cursor to the next feature.
 *
 * Skips blank lines and `#` directives, and stops at the `##FASTA` section that
 * may end a GFF3 file.
 *
 * @param r {Reader} an open cursor
 * @return {Next} `Record` with the feature, or `Done` at end of stream
 * @throws {Error} on a record with fewer than nine columns
 */
export func next(r as Reader) {
    while (source.hasNext($r.src)) {
        def line as string init source.nextLine($r.src);
        if (strings.startsWith($line, "##FASTA")) {
            return Next.Done;
        }
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
 * Parse one GFF3 or GTF record.
 *
 * @param line {string} one nine-column record, newline already stripped
 * @return {Feature} the parsed feature, coordinates left 1-based
 * @throws {Error} when fewer than nine columns are present
 */
export func parse(line as string) {
    try {
        return parseRow($line);
    } catch (e) {
        return core.reraise("gff", $line, $e);
    }
}

/**
 * Parse a GFF/GTF record, assuming the caller has wrapped it for error
 * reporting.
 *
 * @param line {string} one nine-column record, newline already stripped
 * @return {Feature} the parsed feature, coordinates left 1-based
 * @throws {Error} when fewer than nine columns are present
 */
func parseRow(line as string) {
    def f as list of string init core.splitTab($line);
    if (len($f) < 9) {
        def got as string init convert.toString(len($f));
        core.fail("gff", "GFF/GTF record needs 9 columns, got " + $got + ": " + $line);
    }
    def score as float init 0.0;
    if ($f[5] != MISSING and len($f[5]) > 0) {
        $score = convert.toFloat($f[5]);
    }
    def phase as int init -1;
    if ($f[7] != MISSING and len($f[7]) > 0) {
        $phase = convert.toInt($f[7]);
    }
    return Feature{
        seqid: $f[0],
        origin: $f[1],
        kind: $f[2],
        start: convert.toInt($f[3]),
        end: convert.toInt($f[4]),
        score: $score,
        strand: $f[6],
        phase: $phase,
        attributes: parseAttributes($f[8])
    };
}

/**
 * Parse a ninth-column attribute string in either GFF3 or GTF style.
 *
 * The style is detected from the first pair: `key=value` is GFF3, `key "value"`
 * is GTF. GFF3 values are percent-decoded. A column that is empty or the
 * undefined-field marker `.` yields no attributes rather than one oddly-named
 * attribute.
 *
 * @param text {string} the raw ninth column
 * @return {map of string to string} attribute keys to values, in file order
 */
export func parseAttributes(text as string) {
    def out as map of string to string init {};
    if (len(strings.trim($text)) == 0 or strings.trim($text) == MISSING) {
        return $out;
    }
    if (styleOf($text) == STYLE_GTF) {
        return parseGtfAttributes($text);
    }
    for (def pair in strings.split($text, ";")) {
        def trimmed as string init strings.trim($pair);
        if (len($trimmed) == 0) {
            continue;
        }
        def cut as int init strings.indexOf($trimmed, "=");
        if ($cut < 0) {
            $out[$trimmed] = "";
            continue;
        }
        def key as string init strings.substring($trimmed, 0, $cut);
        $out[$key] = percentDecode(strings.substring($trimmed, $cut + 1));
    }
    return $out;
}

/**
 * Which attribute style a ninth column uses.
 *
 * @param text {string} the raw ninth column
 * @return {string} `STYLE_GFF3` or `STYLE_GTF`
 */
export func styleOf(text as string) {
    def first as string init strings.trim(strings.split($text, ";")[0]);
    if (len($first) == 0) {
        return STYLE_GFF3;
    }
    def eq as int init strings.indexOf($first, "=");
    def quote as int init strings.indexOf($first, "\"");
    if ($eq >= 0 and ($quote < 0 or $eq < $quote)) {
        return STYLE_GFF3;
    }
    if ($quote >= 0) {
        return STYLE_GTF;
    }
    # No `=` and no quote: a `key value` pair is GTF with an unquoted value,
    # which GENCODE uses for numeric attributes such as `level 2`.
    if ($eq < 0 and strings.indexOf($first, " ") > 0) {
        return STYLE_GTF;
    }
    return STYLE_GFF3;
}

/**
 * Read an attribute, or a fallback when it is absent.
 *
 * @param f {Feature} the feature
 * @param key {string} attribute name
 * @param fallback {string} value to return when the key is missing
 * @return {string}
 */
export func attribute(f as Feature, key as string, fallback as string) {
    if (not maps.has($f.attributes, $key)) {
        return $fallback;
    }
    return $f.attributes[$key];
}

/**
 * The feature's identifier: GFF3 `ID`, else GTF `transcript_id`, else `gene_id`.
 *
 * @param f {Feature} the feature
 * @return {string} the identifier, "" when none is present
 */
export func identifier(f as Feature) {
    def gff as string init attribute($f, "ID", "");
    if (len($gff) > 0) {
        return $gff;
    }
    def tx as string init attribute($f, "transcript_id", "");
    if (len($tx) > 0) {
        return $tx;
    }
    return attribute($f, "gene_id", "");
}

/**
 * The feature's display name: GFF3 `Name`, else `gene_name`, else `identifier`.
 *
 * @param f {Feature} the feature
 * @return {string} the name, "" when none is present
 */
export func name(f as Feature) {
    def n as string init attribute($f, "Name", "");
    if (len($n) > 0) {
        return $n;
    }
    def gene as string init attribute($f, "gene_name", "");
    if (len($gene) > 0) {
        return $gene;
    }
    return identifier($f);
}

/**
 * Convert a feature to a half-open interval for the interval engine.
 *
 * This is where the 1-based inclusive GFF coordinate becomes the 0-based
 * half-open BED coordinate: `start - 1` and `end` unchanged.
 *
 * @param f {Feature} the feature
 * @return {Interval} the equivalent 0-based half-open interval, named by `name`
 */
export func toInterval(f as Feature) {
    return intervals.named($f.seqid, $f.start - 1, $f.end, name($f), $f.score, $f.strand);
}

/**
 * Build a feature from a half-open interval, converting back to 1-based.
 *
 * @param iv {Interval} a 0-based half-open interval
 * @param origin {string} value for the `source` column
 * @param kind {string} value for the `type` column
 * @return {Feature} the feature, with `ID` set from the interval's name
 */
export func fromInterval(iv as intervals.Interval, origin as string, kind as string) {
    def attrs as map of string to string init {};
    if (len($iv.name) > 0) {
        $attrs["ID"] = $iv.name;
    }
    return Feature{
        seqid: $iv.chrom,
        origin: $origin,
        kind: $kind,
        start: $iv.start + 1,
        end: $iv.end,
        score: $iv.score,
        strand: $iv.strand,
        phase: -1,
        attributes: $attrs
    };
}

/**
 * Render a feature as a nine-column record.
 *
 * @param f {Feature} the feature
 * @param style {string} `STYLE_GFF3` or `STYLE_GTF`
 * @return {string} one tab-delimited record, no trailing newline
 */
export func format(f as Feature, style as string) {
    def phase as string init MISSING;
    if ($f.phase >= 0) {
        $phase = convert.toString($f.phase);
    }
    def score as string init MISSING;
    if ($f.score != 0.0) {
        $score = convert.toString($f.score);
    }
    def columns as list of string init [
        $f.seqid,
        $f.origin,
        $f.kind,
        core.padNum($f.start, 1),
        core.padNum($f.end, 1),
        $score,
        $f.strand,
        $phase,
        formatAttributes($f.attributes, $style)
    ];
    return strings.join($columns, core.TAB);
}

/**
 * Render an attribute map in GFF3 or GTF style.
 *
 * @param attrs {map of string to string} attributes in output order
 * @param style {string} `STYLE_GFF3` or `STYLE_GTF`
 * @return {string} the ninth column, or "." when there are no attributes
 */
export func formatAttributes(attrs as map of string to string, style as string) {
    def parts as list of string init [];
    for (def key in $attrs) {
        if ($style == STYLE_GTF) {
            $parts[] = $key + " \"" + $attrs[$key] + "\"";
        } else {
            $parts[] = $key + "=" + percentEncode($attrs[$key]);
        }
    }
    if (len($parts) == 0) {
        # GFF3 replaces an undefined field with ".", and an empty ninth column
        # is not valid in either dialect.
        return MISSING;
    }
    if ($style == STYLE_GTF) {
        return strings.join($parts, "; ") + ";";
    }
    return strings.join($parts, ";");
}

/**
 * Read a whole annotation file into memory.
 *
 * @param path {string} GFF3/GTF file, optionally gzipped
 * @return {list of Feature} every feature, in file order
 * @throws {Error} when the file is missing or a record is malformed
 */
export func readAll(path as string) {
    def r as Reader init open($path);
    errdefer quietClose($r);
    def out as list of Feature init [];
    # Reads lines directly rather than driving the `Next` cursor, which would
    # copy every feature into an enum payload and back out of it again.
    while (source.hasNext($r.src)) {
        def line as string init source.nextLine($r.src);
        if (strings.startsWith($line, "##FASTA")) {
            break;
        }
        if (not core.isSkippable($line, "#")) {
            $out[] = parse($line);
        }
    }
    close($r);
    return $out;
}

/**
 * Read an annotation file directly as half-open intervals.
 *
 * @param path {string} GFF3/GTF file, optionally gzipped
 * @param kind {string} keep only this feature type, or "" for all types
 * @return {list of Interval} converted features, in file order
 */
export func readIntervals(path as string, kind as string) {
    def out as list of intervals.Interval init [];
    for (def f in readAll($path)) {
        if (len($kind) == 0 or $f.kind == $kind) {
            $out[] = toInterval($f);
        }
    }
    return $out;
}

/**
 * Parse GTF-style `key "value";` attributes.
 *
 * @param text {string} the raw ninth column
 * @return {map of string to string}
 */
func parseGtfAttributes(text as string) {
    def out as map of string to string init {};
    for (def pair in strings.split($text, ";")) {
        def trimmed as string init strings.trim($pair);
        if (len($trimmed) == 0) {
            continue;
        }
        def cut as int init strings.indexOf($trimmed, " ");
        if ($cut < 0) {
            $out[$trimmed] = "";
            continue;
        }
        def key as string init strings.substring($trimmed, 0, $cut);
        $out[$key] = unquote(strings.trim(strings.substring($trimmed, $cut + 1)));
    }
    return $out;
}

/**
 * Strip surrounding double quotes from a GTF attribute value.
 *
 * @param s {string} the raw value
 * @return {string} the value without its quotes
 */
func unquote(s as string) {
    if (len($s) >= 2 and strings.startsWith($s, "\"") and strings.endsWith($s, "\"")) {
        return strings.substring($s, 1, len($s) - 1);
    }
    return $s;
}

/**
 * Percent-decode a GFF3 attribute value.
 *
 * @param s {string} the raw value
 * @return {string} the decoded value, or the input unchanged when it is not
 *   valid percent-encoding
 */
func percentDecode(s as string) {
    if (strings.indexOf($s, "%") < 0) {
        return $s;
    }
    try {
        return convert.stringFromBytes(encoding.fromText($s, "uri-percent"), "utf-8");
    } catch (e) {
        return $s;
    }
}

/**
 * Percent-encode the characters GFF3 reserves in an attribute value.
 *
 * Only the reserved set is escaped (`%` first, so it is not double-escaped):
 * tab, newline, carriage return and `%` per the GFF3 spec's global rule, plus
 * `;`, `=`, `&` and `,` which are reserved in column 9. The spec says "no other
 * characters may be encoded", so spaces are deliberately left alone - GFF3
 * allows them unescaped in attribute values, and encoding them would be wrong.
 *
 * @param s {string} the value to escape
 * @return {string} the escaped value
 */
func percentEncode(s as string) {
    def out as string init strings.replace($s, "%", "%25");
    $out = strings.replace($out, ";", "%3B");
    $out = strings.replace($out, "=", "%3D");
    $out = strings.replace($out, "&", "%26");
    $out = strings.replace($out, ",", "%2C");
    $out = strings.replace($out, core.TAB, "%09");
    $out = strings.replace($out, "\r", "%0D");
    return strings.replace($out, "\n", "%0A");
}

/**
 * Close a cursor on a path that is already failing.
 *
 * Cleanup that runs while an error is propagating must not raise one of its own.
 *
 * @param r {Reader} the cursor to close
 * @return {int} the close result, or -1 when closing itself failed
 */
func quietClose(r as Reader) {
    try {
        return close($r);
    } catch (e) {
        return -1;
    }
}
