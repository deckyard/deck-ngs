# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * White-box tests for the interval engine.
 * Run with: jennifer test intervals_test.j
 */
use testing;

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
