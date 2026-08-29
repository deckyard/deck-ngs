# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0
/**
 * White-box tests for the FASTQ reader.
 * Run with: jennifer test fastq_test.j
 */
use testing;
use fs;

func temp(body as string, suffix as string) {
    def p as string init fs.makeTempFile("", "ngsfastq", $suffix);
    fs.writeString($p, $body);
    return $p;
}

func testParseAndAccessors() {
    def rec as Read init parse("@r1 extra info", "ACGTN", "+", "IIII#");
    testing.assertEqual($rec.id, "r1");
    testing.assertEqual($rec.desc, "extra info");
    testing.assertEqual(length($rec), 5);
    testing.assertEqual(ambiguousCount($rec), 1);
    testing.assertEqual(gcPercent($rec), 40.0);
    testing.assertEqual(minQuality($rec, PHRED33), 2);
    testing.assertEqual(headerOf($rec), "r1 extra info");
}

func testQualityDecoding() {
    def rec as Read init parse("@r", "AC", "+", "I!");
    def qs as list of int init qualities($rec, PHRED33);
    testing.assertEqual($qs[0], 40);
    testing.assertEqual($qs[1], 0);
    testing.assertEqual(meanQuality($rec, PHRED33), 20.0);
}

func testTruncatedRecordFails() {
    testing.assertThrows("truncatedFastq", "ngs.fastq");
}

func truncatedFastq() {
    def p as string init temp("@a\nACGT\n+\nIIII\n@b\nACGT\n", ".fastq");
    def r as Reader init open($p);
    def n as int init 0;
    while (hasNext($r)) {
        $n = $n + length(next($r));
    }
    return $n;
}

func testLengthMismatchFails() {
    testing.assertThrows("badFastq", "ngs.fastq");
}

func badFastq() {
    return parse("@r", "ACGT", "+", "II");
}

func testStreamAndFormat() {
    def p as string init temp("@a\nACGT\n+\nIIII\n@b desc\nTTTT\n+\n####\n", ".fastq");
    def r as Reader init open($p);
    def got as list of string init [];
    while (hasNext($r)) {
        def rec as Read init next($r);
        $got[] = $rec.id;
    }
    closeStrict($r);
    fs.remove($p);
    testing.assertEqual(len($got), 2);
    testing.assertEqual($got[1], "b");
}

func testRoundTripFormat() {
    def rec as Read init parse("@r1 note", "ACGT", "+", "IIII");
    testing.assertEqual(format($rec), "@r1 note\nACGT\n+\nIIII\n");
    testing.assertEqual(toFasta($rec), ">r1 note\nACGT\n");
}

func testSubseqAndRevComp() {
    def rec as Read init parse("@r", "ACGTA", "+", "12345");
    def cut as Read init subseq($rec, 1, 4);
    testing.assertEqual($cut.seq, "CGT");
    testing.assertEqual($cut.qual, "234");
    def rc as Read init reverseComplement($rec);
    testing.assertEqual($rc.seq, "TACGT");
    testing.assertEqual($rc.qual, "54321");
}
