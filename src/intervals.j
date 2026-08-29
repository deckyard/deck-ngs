# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Interval arithmetic over genomic features - the bedtools-lite layer.
 *
 * Coordinates are **BED convention: 0-based, half-open** `[start, end)`. GFF and
 * VCF are 1-based inclusive, so `gff.toInterval` and `vcf.toInterval` do the
 * conversion for you; never mix the two by hand.
 *
 * Sorting is lexicographic by chromosome then numeric by start, matching
 * `sort -k1,1 -k2,2n` and bedtools' default order - so `chr10` sorts before
 * `chr2`.
 *
 * The set operations (`intersect`, `subtract`, `coverage`, `countOverlaps`) run
 * a sorted sweep rather than a nested scan, so they are linear in the input plus
 * the number of reported overlaps.
 *
 * @module intervals
 * @since 0.1.0
 */
use lists;
import "./core.j" as core;

/** Width used when zero-padding coordinates into sort keys. */
export def const KEY_WIDTH as int init 12;

/**
 * One past the largest coordinate an interval may carry.
 *
 * Sorting compares zero-padded coordinates as text, so a value too wide for
 * `KEY_WIDTH` would pad to nothing and sort by its first digit - putting
 * 1000000000000 before 999999999999. Rejecting it at construction turns a silent
 * misordering into an error at the point the bad value entered.
 *
 * The bound is 10^12 bases, about 750 times the largest chromosome ever
 * sequenced, so no real coordinate approaches it.
 */
export def const MAX_COORDINATE as int init 1000000000000;

/**
 * A half-open genomic interval.
 *
 * @field chrom {string} reference sequence name
 * @field start {int} 0-based inclusive start
 * @field end {int} 0-based exclusive end
 * @field name {string} feature name, "" when unnamed
 * @field score {float} feature score, 0.0 when absent
 * @field strand {string} "+", "-" or "." when unstranded
 */
export def struct Interval {
    chrom as string,
    start as int,
    end as int,
    name as string,
    score as float,
    strand as string
};

/**
 * Build an unnamed, unstranded interval.
 *
 * @param chrom {string} reference sequence name
 * @param start {int} 0-based inclusive start
 * @param end {int} 0-based exclusive end
 * @return {Interval}
 * @throws {Error} when the coordinates are negative or inverted
 */
export func make(chrom as string, start as int, end as int) {
    return named($chrom, $start, $end, "", 0.0, ".");
}

/**
 * Build a fully specified interval.
 *
 * @param chrom {string} reference sequence name
 * @param start {int} 0-based inclusive start
 * @param end {int} 0-based exclusive end
 * @param name {string} feature name
 * @param score {float} feature score
 * @param strand {string} "+", "-" or "."
 * @return {Interval}
 * @throws {Error} when the coordinates are negative, inverted, or too large to
 *   sort correctly (see `MAX_COORDINATE`)
 */
export func named(
    chrom as string,
    start as int,
    end as int,
    name as string,
    score as float,
    strand as string) {
    if ($start < 0 or $end < 0) {
        core.fail("intervals", "negative coordinate in " + $chrom);
    }
    if ($end < $start) {
        core.fail(
            "intervals",
            "inverted interval " + $chrom + ":" + core.padNum($start, 1) + "-" +
                core.padNum($end, 1));
    }
    if ($end >= MAX_COORDINATE) {
        core.fail(
            "intervals",
            "coordinate too large to sort correctly (limit " +
                core.padNum(MAX_COORDINATE, 1) + "): " + $chrom + ":" + core.padNum($end, 1));
    }
    return Interval{
        chrom: $chrom,
        start: $start,
        end: $end,
        name: $name,
        score: $score,
        strand: $strand
    };
}

/**
 * Width of an interval in bases.
 *
 * @param iv {Interval}
 * @return {int}
 */
export func length(iv as Interval) {
    return $iv.end - $iv.start;
}

/**
 * Whether two intervals share at least one base.
 *
 * Book-ended intervals (`a.end == b.start`) do not overlap.
 *
 * @param a {Interval}
 * @param b {Interval}
 * @return {bool}
 */
export func overlaps(a as Interval, b as Interval) {
    if ($a.chrom != $b.chrom) {
        return false;
    }
    return $a.start < $b.end and $b.start < $a.end;
}

/**
 * Number of bases two intervals share.
 *
 * @param a {Interval}
 * @param b {Interval}
 * @return {int} 0 when they do not overlap
 */
export func overlapLength(a as Interval, b as Interval) {
    if (not overlaps($a, $b)) {
        return 0;
    }
    return lowerOf($a.end, $b.end) - upperOf($a.start, $b.start);
}

/**
 * Whether `outer` fully contains `inner`.
 *
 * @param outer {Interval}
 * @param inner {Interval}
 * @return {bool}
 */
export func contains(outer as Interval, inner as Interval) {
    if ($outer.chrom != $inner.chrom) {
        return false;
    }
    return $outer.start <= $inner.start and $inner.end <= $outer.end;
}

/**
 * The shared region of two overlapping intervals.
 *
 * @param a {Interval}
 * @param b {Interval}
 * @return {Interval} the intersection, named after `a`
 * @throws {Error} when the two do not overlap
 */
export func intersection(a as Interval, b as Interval) {
    if (not overlaps($a, $b)) {
        core.fail("intervals", "cannot intersect non-overlapping intervals");
    }
    return named(
        $a.chrom,
        upperOf($a.start, $b.start),
        lowerOf($a.end, $b.end),
        $a.name,
        $a.score,
        $a.strand);
}

/**
 * Sort key for one interval: chromosome, then zero-padded start and end.
 *
 * Exposed because `lists.sortBy` needs it as a `func` value; you rarely call it
 * directly.
 *
 * @param iv {Interval}
 * @return {string} a lexicographically sortable key
 */
export func sortKey(iv as Interval) {
    return $iv.chrom + core.TAB + core.padNum($iv.start, KEY_WIDTH) + core.TAB +
        core.padNum($iv.end, KEY_WIDTH);
}

/**
 * Sort intervals by chromosome then start then end.
 *
 * Input that is already in order is returned as-is rather than re-sorted, which
 * is the common case for BED and GFF files coming out of a pipeline. Value
 * semantics mean the caller still binds its own copy either way, so this is not
 * observable beyond being faster.
 *
 * @param ivs {list of Interval} intervals to sort
 * @return {list of Interval} the intervals in sort order
 */
export func sort(ivs as list of Interval) {
    if (isSorted($ivs)) {
        return $ivs;
    }
    return lists.sortBy($ivs, sortKey);
}

/**
 * Whether a list is already in sort order.
 *
 * Compares fields directly, so it builds no keys at all: a linear pass of
 * integer and string comparisons against `n` key strings plus a sort. Real BED
 * and GFF files usually arrive coordinate-sorted, and every set operation sorts
 * its inputs, so this turns the common case from `O(n log n)` plus `n`
 * allocations into one scan. Unsorted input costs almost nothing either, since
 * the scan stops at the first pair out of order.
 *
 * @param ivs {list of Interval} the list to inspect
 * @return {bool} true when no element precedes the one before it
 */
func isSorted(ivs as list of Interval) {
    if (len($ivs) < 2) {
        return true;
    }
    # Indexes the list directly rather than binding each element: binding an
    # element of a `list of Interval` to a local deep-copies the struct, which
    # over a large annotation is the dominant cost of the scan.
    for (def i in 1..len($ivs)) {
        if ($ivs[$i].chrom < $ivs[$i - 1].chrom) {
            return false;
        }
        if ($ivs[$i].chrom == $ivs[$i - 1].chrom) {
            if ($ivs[$i].start < $ivs[$i - 1].start) {
                return false;
            }
            if ($ivs[$i].start == $ivs[$i - 1].start and $ivs[$i].end < $ivs[$i - 1].end) {
                return false;
            }
        }
    }
    return true;
}

/**
 * Merge overlapping and nearby intervals, per chromosome.
 *
 * @param ivs {list of Interval} intervals in any order
 * @param maxGap {int} merge features separated by at most this many bases; 0
 *   merges overlapping and book-ended features, matching `bedtools merge -d 0`
 * @return {list of Interval} sorted, disjoint intervals with names cleared
 */
export func merge(ivs as list of Interval, maxGap as int) {
    def ordered as list of Interval init sort($ivs);
    def out as list of Interval init [];
    if (len($ordered) == 0) {
        return $out;
    }
    def cur as Interval init $ordered[0];
    for (def i in 1..len($ordered)) {
        # Reads fields off the indexed element instead of binding it: a binding
        # would deep-copy the struct once per element, which is what made a
        # merge of a large annotation allocate more than the file itself.
        if ($ordered[$i].chrom == $cur.chrom and $ordered[$i].start <= $cur.end + $maxGap) {
            # Extending the run is a single field write; rebuilding it through
            # `make` cost a call and a re-validation per element.
            if ($ordered[$i].end > $cur.end) {
                $cur.end = $ordered[$i].end;
            }
        } else {
            $out[] = bare($cur);
            $cur = $ordered[$i];
        }
    }
    $out[] = bare($cur);
    return $out;
}

/**
 * Total number of bases covered, counting overlaps once.
 *
 * @param ivs {list of Interval} intervals in any order
 * @return {int} bases covered by the union of the input
 */
export func totalLength(ivs as list of Interval) {
    def merged as list of Interval init merge($ivs, 0);
    def total as int init 0;
    for (def i in 0..len($merged)) {
        $total = $total + ($merged[$i].end - $merged[$i].start);
    }
    return $total;
}

/**
 * Pairwise intersection regions between two interval sets.
 *
 * Equivalent to `bedtools intersect -a A -b B`: one output interval per
 * overlapping pair, covering only the shared bases and carrying `a`'s name.
 *
 * @param a {list of Interval} query intervals
 * @param b {list of Interval} target intervals
 * @return {list of Interval} shared regions, sorted
 */
export func intersect(a as list of Interval, b as list of Interval) {
    def qs as list of Interval init sort($a);
    def ts as list of Interval init sort($b);
    def out as list of Interval init [];
    def cursor as int init 0;
    for (def q in $qs) {
        $cursor = advanceTo($ts, $cursor, $q);
        def k as int init $cursor;
        while ($k < len($ts) and $ts[$k].chrom == $q.chrom and $ts[$k].start < $q.end) {
            if ($q.start < $ts[$k].end) {
                $out[] = named(
                    $q.chrom,
                    upperOf($q.start, $ts[$k].start),
                    lowerOf($q.end, $ts[$k].end),
                    $q.name,
                    $q.score,
                    $q.strand);
            }
            $k = $k + 1;
        }
    }
    return $out;
}

/**
 * Count, for each query interval, how many target intervals it overlaps.
 *
 * Equivalent to `bedtools intersect -c`.
 *
 * @param a {list of Interval} query intervals
 * @param b {list of Interval} target intervals
 * @return {list of int} one count per interval of `sort(a)`, in that order
 */
export func countOverlaps(a as list of Interval, b as list of Interval) {
    def qs as list of Interval init sort($a);
    def ts as list of Interval init sort($b);
    def out as list of int init [];
    def cursor as int init 0;
    for (def q in $qs) {
        $cursor = advanceTo($ts, $cursor, $q);
        def k as int init $cursor;
        def hits as int init 0;
        while ($k < len($ts) and $ts[$k].chrom == $q.chrom and $ts[$k].start < $q.end) {
            if ($q.start < $ts[$k].end) {
                $hits = $hits + 1;
            }
            $k = $k + 1;
        }
        $out[] = $hits;
    }
    return $out;
}

/**
 * Bases of each query interval covered by the target set.
 *
 * Overlapping targets are merged first, so a base covered twice counts once.
 * Equivalent to the covered-bases column of `bedtools coverage`.
 *
 * @param a {list of Interval} query intervals
 * @param b {list of Interval} target intervals
 * @return {list of int} covered base counts, aligned with `sort(a)`
 */
export func coverage(a as list of Interval, b as list of Interval) {
    def qs as list of Interval init sort($a);
    def ts as list of Interval init merge($b, 0);
    def out as list of int init [];
    def cursor as int init 0;
    for (def q in $qs) {
        $cursor = advanceTo($ts, $cursor, $q);
        def k as int init $cursor;
        def covered as int init 0;
        while ($k < len($ts) and $ts[$k].chrom == $q.chrom and $ts[$k].start < $q.end) {
            if ($q.start < $ts[$k].end) {
                $covered = $covered +
                    lowerOf($q.end, $ts[$k].end) - upperOf($q.start, $ts[$k].start);
            }
            $k = $k + 1;
        }
        $out[] = $covered;
    }
    return $out;
}

/**
 * Remove from each query interval every base covered by the target set.
 *
 * Equivalent to `bedtools subtract`. A query fully covered disappears; a query
 * split in the middle yields two pieces.
 *
 * @param a {list of Interval} query intervals
 * @param b {list of Interval} intervals to remove
 * @return {list of Interval} the remaining pieces, sorted, names preserved
 */
export func subtract(a as list of Interval, b as list of Interval) {
    def qs as list of Interval init sort($a);
    def ts as list of Interval init merge($b, 0);
    def out as list of Interval init [];
    def cursor as int init 0;
    for (def q in $qs) {
        $cursor = advanceTo($ts, $cursor, $q);
        def pieces as list of Interval init clipPieces($q, $ts, $cursor);
        for (def piece in $pieces) {
            $out[] = $piece;
        }
    }
    return $out;
}

/**
 * The parts of one query interval left uncovered by a merged target set.
 *
 * @param q {Interval} the query interval
 * @param ts {list of Interval} sorted, disjoint targets
 * @param from {int} sweep cursor into `ts`
 * @return {list of Interval} the surviving pieces of `q`, in order
 */
func clipPieces(q as Interval, ts as list of Interval, from as int) {
    def pieces as list of Interval init [];
    def at as int init $q.start;
    def k as int init $from;
    while ($k < len($ts) and $ts[$k].chrom == $q.chrom and $ts[$k].start < $q.end) {
        def edge as int init $ts[$k].start;
        def stop as int init $ts[$k].end;
        $k = $k + 1;
        if ($stop <= $at or $edge >= $q.end) {
            continue;
        }
        if ($edge > $at) {
            $pieces[] = named($q.chrom, $at, $edge, $q.name, $q.score, $q.strand);
        }
        $at = $stop;
    }
    if ($at < $q.end) {
        $pieces[] = named($q.chrom, $at, $q.end, $q.name, $q.score, $q.strand);
    }
    return $pieces;
}

/**
 * The gaps between intervals, across whole chromosomes.
 *
 * Equivalent to `bedtools complement -g`.
 *
 * @param ivs {list of Interval} covered intervals
 * @param sizes {map of string to int} chromosome name to length
 * @return {list of Interval} uncovered regions, sorted
 */
export func complement(ivs as list of Interval, sizes as map of string to int) {
    def merged as list of Interval init merge($ivs, 0);
    def out as list of Interval init [];
    for (def chrom in $sizes) {
        def gaps as list of Interval init gapsIn($merged, $chrom, $sizes[$chrom]);
        for (def gap in $gaps) {
            $out[] = $gap;
        }
    }
    return sort($out);
}

/**
 * Uncovered stretches of one chromosome, given merged covered intervals.
 *
 * @param merged {list of Interval} sorted, disjoint covered intervals
 * @param chrom {string} the chromosome to walk
 * @param limit {int} the chromosome's length
 * @return {list of Interval} the gaps, in order
 */
func gapsIn(merged as list of Interval, chrom as string, limit as int) {
    def out as list of Interval init [];
    def at as int init 0;
    for (def i in 0..len($merged)) {
        if ($merged[$i].chrom != $chrom) {
            continue;
        }
        if ($merged[$i].start > $at) {
            $out[] = make($chrom, $at, lowerOf($merged[$i].start, $limit));
        }
        if ($merged[$i].end > $at) {
            $at = $merged[$i].end;
        }
    }
    if ($at < $limit) {
        $out[] = make($chrom, $at, $limit);
    }
    return $out;
}

/**
 * Grow an interval on both sides, clamped at zero and an optional chromosome end.
 *
 * Equivalent to `bedtools slop`.
 *
 * @param iv {Interval} the interval to grow
 * @param left {int} bases to add before the start
 * @param right {int} bases to add after the end
 * @param limit {int} chromosome length to clamp against; 0 or less means no clamp
 * @return {Interval} the grown interval
 */
export func slop(iv as Interval, left as int, right as int, limit as int) {
    def start as int init $iv.start - $left;
    if ($start < 0) {
        $start = 0;
    }
    def end as int init $iv.end + $right;
    if ($limit > 0 and $end > $limit) {
        $end = $limit;
    }
    return named($iv.chrom, $start, $end, $iv.name, $iv.score, $iv.strand);
}

/**
 * Group intervals by chromosome.
 *
 * The input is sorted first, so each group is itself sorted by start and the
 * grouping costs one pass rather than a copy-back per element.
 *
 * @param ivs {list of Interval} intervals in any order
 * @return {map of string to list of Interval} sorted intervals per chromosome
 */
export func byChromosome(ivs as list of Interval) {
    def ordered as list of Interval init sort($ivs);
    def out as map of string to list of Interval init {};
    def current as string init "";
    def bucket as list of Interval init [];
    for (def iv in $ordered) {
        if ($iv.chrom != $current) {
            if (len($bucket) > 0) {
                $out[$current] = $bucket;
            }
            $current = $iv.chrom;
            $bucket = [];
        }
        $bucket[] = $iv;
    }
    if (len($bucket) > 0) {
        $out[$current] = $bucket;
    }
    return $out;
}

/**
 * Render an interval as a three-column BED line.
 *
 * @param iv {Interval}
 * @return {string} tab-delimited `chrom start end`, no trailing newline
 */
export func toBed3(iv as Interval) {
    return $iv.chrom + core.TAB + core.padNum($iv.start, 1) + core.TAB + core.padNum($iv.end, 1);
}

/**
 * A copy of an interval with its name, score and strand cleared.
 *
 * A merged run is not any one of its inputs, so carrying one input's name onto
 * it would be a lie. Built as a literal rather than through `make`, which would
 * re-validate coordinates that came from already-validated intervals.
 *
 * @param iv {Interval} the interval to strip
 * @return {Interval} the same span, unnamed and unstranded
 */
func bare(iv as Interval) {
    return Interval{
        chrom: $iv.chrom,
        start: $iv.start,
        end: $iv.end,
        name: "",
        score: 0.0,
        strand: "."
    };
}

/**
 * Advance a sweep cursor past every target that ends at or before the query start.
 *
 * Safe because both lists are sorted by start: a target ending before this
 * query's start cannot overlap any later query either.
 *
 * @param ts {list of Interval} sorted targets
 * @param from {int} current cursor position
 * @param q {Interval} the query being processed
 * @return {int} the new cursor position
 */
func advanceTo(ts as list of Interval, from as int, q as Interval) {
    def i as int init $from;
    while ($i < len($ts)) {
        if ($ts[$i].chrom != $q.chrom) {
            if ($ts[$i].chrom > $q.chrom) {
                return $i;
            }
        } elseif ($ts[$i].end > $q.start) {
            return $i;
        }
        $i = $i + 1;
    }
    return $i;
}

/**
 * The smaller of two ints.
 *
 * @param a {int}
 * @param b {int}
 * @return {int}
 */
func lowerOf(a as int, b as int) {
    if ($a < $b) {
        return $a;
    }
    return $b;
}

/**
 * The larger of two ints.
 *
 * @param a {int}
 * @param b {int}
 * @return {int}
 */
func upperOf(a as int, b as int) {
    if ($a > $b) {
        return $a;
    }
    return $b;
}
