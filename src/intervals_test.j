# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

/**
 * White-box tests for the interval engine.
 * Run with: jennifer test intervals_test.j
 */
use testing;
use math;
use time;
use convert;

func iv(c as string, s as int, e as int) {
    return make($c, $s, $e);
}

func testLengthAndOverlap() {
    testing.assertEqual(length(iv("chr1", 10, 20)), 10);
    testing.assertEqual(overlaps(iv("chr1", 10, 20), iv("chr1", 15, 25)), true);
    testing.assertEqual(overlaps(iv("chr1", 10, 20), iv("chr1", 20, 25)), false);
    testing.assertEqual(overlaps(iv("chr1", 10, 20), iv("chr2", 15, 25)), false);
    testing.assertEqual(overlapLength(iv("chr1", 10, 20), iv("chr1", 15, 25)), 5);
    testing.assertEqual(overlapLength(iv("chr1", 10, 20), iv("chr1", 30, 40)), 0);
    testing.assertEqual(contains(iv("chr1", 10, 30), iv("chr1", 15, 25)), true);
    testing.assertEqual(contains(iv("chr1", 10, 30), iv("chr1", 5, 25)), false);
}

func testSortIsLexicographicThenNumeric() {
    def xs as list of Interval init [
        iv("chr2", 5, 6),
        iv("chr10", 1, 2),
        iv("chr1", 100, 200),
        iv("chr1", 20, 30)
    ];
    def got as list of Interval init sort($xs);
    testing.assertEqual($got[0].chrom, "chr1");
    testing.assertEqual($got[0].start, 20);
    testing.assertEqual($got[1].start, 100);
    testing.assertEqual($got[2].chrom, "chr10");
    testing.assertEqual($got[3].chrom, "chr2");
}

func testMergeOverlappingAndBookended() {
    def xs as list of Interval init [
        iv("chr1", 10, 20),
        iv("chr1", 15, 25),
        iv("chr1", 25, 30),
        iv("chr1", 40, 50),
        iv("chr2", 1, 5)
    ];
    def got as list of Interval init merge($xs, 0);
    testing.assertEqual(len($got), 3);
    testing.assertEqual($got[0].start, 10);
    testing.assertEqual($got[0].end, 30);
    testing.assertEqual($got[1].start, 40);
    testing.assertEqual($got[2].chrom, "chr2");
}

func testMergeWithGap() {
    def xs as list of Interval init [iv("chr1", 10, 20), iv("chr1", 25, 30)];
    testing.assertEqual(len(merge($xs, 0)), 2);
    testing.assertEqual(len(merge($xs, 5)), 1);
    testing.assertEqual(merge($xs, 5)[0].end, 30);
}

func testTotalLengthCountsUnionOnce() {
    def xs as list of Interval init [iv("chr1", 0, 100), iv("chr1", 50, 150)];
    testing.assertEqual(totalLength($xs), 150);
}

func testIntersectReportsSharedRegions() {
    def a as list of Interval init [iv("chr1", 10, 30), iv("chr1", 100, 200)];
    def b as list of Interval init [
        iv("chr1", 20, 25),
        iv("chr1", 28, 40),
        iv("chr1", 150, 300),
        iv("chr2", 0, 999)
    ];
    def got as list of Interval init intersect($a, $b);
    testing.assertEqual(len($got), 3);
    testing.assertEqual($got[0].start, 20);
    testing.assertEqual($got[0].end, 25);
    testing.assertEqual($got[1].start, 28);
    testing.assertEqual($got[1].end, 30);
    testing.assertEqual($got[2].start, 150);
    testing.assertEqual($got[2].end, 200);
}

func testIntersectHandlesOneQueryManyTargets() {
    def a as list of Interval init [iv("chr1", 0, 1000)];
    def b as list of Interval init [iv("chr1", 10, 20), iv("chr1", 30, 40), iv("chr1", 50, 60)];
    testing.assertEqual(len(intersect($a, $b)), 3);
}

func testCountOverlaps() {
    def a as list of Interval init [iv("chr1", 0, 50), iv("chr1", 60, 70)];
    def b as list of Interval init [iv("chr1", 10, 20), iv("chr1", 30, 40), iv("chr1", 100, 110)];
    def got as list of int init countOverlaps($a, $b);
    testing.assertEqual($got[0], 2);
    testing.assertEqual($got[1], 0);
}

func testCoverageCountsEachBaseOnce() {
    def a as list of Interval init [iv("chr1", 0, 100)];
    def b as list of Interval init [iv("chr1", 10, 30), iv("chr1", 20, 40)];
    testing.assertEqual(coverage($a, $b)[0], 30);
}

func testSubtractSplitsAndRemoves() {
    def a as list of Interval init [iv("chr1", 0, 100), iv("chr1", 200, 210)];
    def b as list of Interval init [iv("chr1", 40, 60), iv("chr1", 200, 210)];
    def got as list of Interval init subtract($a, $b);
    testing.assertEqual(len($got), 2);
    testing.assertEqual($got[0].start, 0);
    testing.assertEqual($got[0].end, 40);
    testing.assertEqual($got[1].start, 60);
    testing.assertEqual($got[1].end, 100);
}

func testSubtractLeavesUntouchedQuery() {
    def a as list of Interval init [iv("chr1", 0, 10)];
    def b as list of Interval init [iv("chr2", 0, 10)];
    def got as list of Interval init subtract($a, $b);
    testing.assertEqual(len($got), 1);
    testing.assertEqual($got[0].end, 10);
}

func testComplementFillsGaps() {
    def xs as list of Interval init [iv("chr1", 10, 20), iv("chr1", 30, 40)];
    def sizes as map of string to int init {"chr1": 100};
    def got as list of Interval init complement($xs, $sizes);
    testing.assertEqual(len($got), 3);
    testing.assertEqual($got[0].start, 0);
    testing.assertEqual($got[0].end, 10);
    testing.assertEqual($got[1].start, 20);
    testing.assertEqual($got[1].end, 30);
    testing.assertEqual($got[2].start, 40);
    testing.assertEqual($got[2].end, 100);
}

func testSlopClamps() {
    def got as Interval init slop(iv("chr1", 5, 10), 20, 20, 25);
    testing.assertEqual($got.start, 0);
    testing.assertEqual($got.end, 25);
    def free as Interval init slop(iv("chr1", 5, 10), 2, 3, 0);
    testing.assertEqual($free.start, 3);
    testing.assertEqual($free.end, 13);
}

func testByChromosome() {
    def xs as list of Interval init [iv("chr2", 5, 6), iv("chr1", 1, 2), iv("chr1", 0, 1)];
    def got as map of string to list of Interval init byChromosome($xs);
    testing.assertEqual(len($got["chr1"]), 2);
    testing.assertEqual($got["chr1"][0].start, 0);
    testing.assertEqual(len($got["chr2"]), 1);
}

func testEmptyInputs() {
    def none as list of Interval init [];
    testing.assertEqual(len(merge($none, 0)), 0);
    testing.assertEqual(len(intersect($none, $none)), 0);
    testing.assertEqual(totalLength($none), 0);
}

func testCoordinateBeyondSortWidthFails() {
    # Sorting pads coordinates to KEY_WIDTH and compares as text, so a value too
    # wide would silently sort before narrower ones. It must be refused instead.
    testing.assertEqual(length(make("chr1", 999999999998, 999999999999)), 1);
    testing.assertThrows("buildOversized", "ngs.intervals");
}

func buildOversized() {
    return make("chr1", 0, MAX_COORDINATE);
}

func testInvertedIntervalFails() {
    testing.assertThrows("buildInverted", "ngs.intervals");
}

func buildInverted() {
    return make("chr1", 50, 10);
}

func testIntersectionKeepsTheSharedBasesAndANames() {
    def a as Interval init named("chr1", 10, 50, "a", 3.0, "+");
    def b as Interval init named("chr1", 30, 80, "b", 9.0, "-");
    def x as Interval init intersection($a, $b);
    testing.assertEqual($x.chrom, "chr1");
    testing.assertEqual($x.start, 30);
    testing.assertEqual($x.end, 50);
    testing.assertEqual($x.name, "a");
    testing.assertEqual($x.score, 3.0);
    testing.assertEqual($x.strand, "+");
    def inner as Interval init intersection($b, iv("chr1", 40, 45));
    testing.assertEqual(length($inner), 5);
    testing.assertEqual($inner.name, "b");
}

func testIntersectionOfBookendedIntervalsFails() {
    testing.assertThrows("intersectBookended", "ngs.intervals");
    testing.assertThrows("intersectAcrossChromosomes", "ngs.intervals");
}

func intersectBookended() {
    return intersection(iv("chr1", 10, 50), iv("chr1", 50, 60));
}

func intersectAcrossChromosomes() {
    return intersection(iv("chr1", 10, 50), iv("chr2", 10, 50));
}

func testSortKeyPadsCoordinatesToKeyWidth() {
    testing.assertEqual(sortKey(iv("chr1", 10, 50)), "chr1\t000000000010\t000000000050");
    # Padding is what makes 9 sort before 10 as text.
    testing.assertEqual(sortKey(iv("chr1", 9, 50)) < sortKey(iv("chr1", 10, 50)), true);
    testing.assertEqual(sortKey(iv("chr1", 10, 9000)) < sortKey(iv("chr1", 10, 10000)), true);
}

func testToBed3() {
    testing.assertEqual(toBed3(named("chrX", 0, 155270560, "x", 1.0, "-")), "chrX\t0\t155270560");
}

func nowMillis() {
    return time.unixMillis(time.now());
}

# Random sets on two chromosomes with long, short and zero-length intervals;
# the expected pairs come from comparing every query with every target.
func randomSet(n as int) {
    def out as list of Interval init [];
    for (def i in 0..$n) {
        def chrom as string init "chr1";
        if (math.randInt(0, 2) == 0) {
            $chrom = "chr2";
        }
        def start as int init math.randInt(0, 1000);
        def width as int init math.randInt(0, 30);
        if (math.randInt(0, 9) == 0) {
            $width = math.randInt(200, 900);
        }
        $out[] = named($chrom, $start, $start + $width, "q" + convert.toString($i), 0.0, "+");
    }
    return $out;
}

# Every target `q` overlaps, as its intersection with `q`, by direct comparison.
func bruteOverlaps(q as Interval, ts as list of Interval) {
    def out as list of Interval init [];
    for (def j in 0..len($ts)) {
        if (overlaps($q, $ts[$j])) {
            $out[] = intersection($q, $ts[$j]);
        }
    }
    return $out;
}

func testIntersectAndCountOverlapsMatchBruteForce() {
    math.randSeed(7);
    for (def round in 0..5) {
        def qs as list of Interval init sort(randomSet(120));
        def ts as list of Interval init sort(randomSet(150));
        def want as list of Interval init [];
        def counts as list of int init [];
        for (def i in 0..len($qs)) {
            def hits as list of Interval init bruteOverlaps($qs[$i], $ts);
            $counts[] = len($hits);
            for (def h in $hits) {
                $want[] = $h;
            }
        }
        testing.assertEqual(intersect($qs, $ts), $want);
        testing.assertEqual(countOverlaps($qs, $ts), $counts);
    }
}

func testOneLongTargetDoesNotMakeIntersectQuadratic() {
    # Each query [10i+5, 10i+6) misses the point target [10i, 10i+1) and lies
    # inside the long one, so every count is 1. Before the active-list sweep the
    # long target pinned the cursor and n = 3000 took about 7 s.
    def n as int init 3000;
    def ts as list of Interval init [iv("chr1", 0, 900000000)];
    def qs as list of Interval init [];
    for (def i in 0..$n) {
        $ts[] = iv("chr1", $i * 10, $i * 10 + 1);
        $qs[] = iv("chr1", $i * 10 + 5, $i * 10 + 6);
    }
    def t0 as int init nowMillis();
    def counts as list of int init countOverlaps($qs, $ts);
    def shared as list of Interval init intersect($qs, $ts);
    def took as int init nowMillis() - $t0;
    testing.assertEqual(len($shared), $n);
    testing.assertEqual($counts[0], 1);
    testing.assertEqual($counts[$n - 1], 1);
    testing.assertEqual($shared[$n - 1].start, ($n - 1) * 10 + 5);
    testing.assertEqual($took < 2000, true);
}

func testComplementClipsIntervalsPastTheChromosomeEnd() {
    def sizes as map of string to int init {"chr1": 100};
    # Both intervals start past the end: the whole chromosome is uncovered.
    def beyond as list of Interval init complement(
        [iv("chr1", 150, 160), iv("chr1", 170, 180)],
        $sizes);
    testing.assertEqual(len($beyond), 1);
    testing.assertEqual($beyond[0].start, 0);
    testing.assertEqual($beyond[0].end, 100);
    # One interval straddles the end: only [0, 50) is left.
    def straddle as list of Interval init complement([iv("chr1", 50, 160)], $sizes);
    testing.assertEqual(len($straddle), 1);
    testing.assertEqual($straddle[0].end, 50);
}

func testComplementReturnsBareChromosomesInSortOrder() {
    def sizes as map of string to int init {"chr2": 50, "chr1": 30, "chr10": 20};
    def got as list of Interval init complement([iv("chr1", 10, 20), iv("chrUn", 0, 5)], $sizes);
    # chr1 has two gaps; chr10 and chr2 have no intervals and come back whole,
    # in sort order (chr10 before chr2). chrUn is not in `sizes`.
    testing.assertEqual(
        $got,
        [iv("chr1", 0, 10), iv("chr1", 20, 30), iv("chr10", 0, 20), iv("chr2", 0, 50)]);
}

func testComplementIsLinearInContigs() {
    # 3000 contigs of 1000 bases, each covered at [10, 20): two gaps per contig.
    # Walking every interval once per contig took about 6 s at this size.
    def sizes as map of string to int init {};
    def ivs as list of Interval init [];
    for (def i in 0..3000) {
        def name as string init "ctg" + convert.toString($i);
        $sizes[$name] = 1000;
        $ivs[] = iv($name, 10, 20);
    }
    def t0 as int init nowMillis();
    def got as list of Interval init complement($ivs, $sizes);
    def took as int init nowMillis() - $t0;
    testing.assertEqual(len($got), 6000);
    testing.assertEqual($took < 2000, true);
}

func testNamedAcceptsOnlyKnownStrands() {
    testing.assertEqual(named("chr1", 0, 1, "", 0.0, "+").strand, "+");
    testing.assertEqual(named("chr1", 0, 1, "", 0.0, "-").strand, "-");
    testing.assertEqual(named("chr1", 0, 1, "", 0.0, ".").strand, ".");
    testing.assertThrows("buildBadStrand", "ngs.intervals");
    testing.assertThrows("buildEmptyStrand", "ngs.intervals");
}

func buildBadStrand() {
    return named("chr1", 0, 1, "", 0.0, "x");
}

func buildEmptyStrand() {
    return named("chr1", 0, 1, "", 0.0, "");
}
