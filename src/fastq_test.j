# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
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

# New fixtures live inside the checkout, under the git-ignored .ci-fixtures/:
# in CI the interpreter runs in a container that mounts only the workspace.
func fixtureDir() {
    fs.mkdirAll(".ci-fixtures");
    return fs.makeTempDir(".ci-fixtures", "ngsfastq");
}

func testWriteRoundTrips() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def p as string init $dir + "/out.fastq";
    def k as source.Sink init source.createPlain($p);
    write($k, parse("@r1 lane 1", "ACGT", "+", "II#I"));
    write($k, parse("@r2", "GG", "+r2", "!!"));
    testing.assertEqual(source.closeSink($k), 0);
    # The separator is always written bare, whatever the input carried.
    testing.assertEqual(fs.readString($p), "@r1 lane 1\nACGT\n+\nII#I\n@r2\nGG\n+\n!!\n");
    def r as Reader init open($p);
    def first as Read init next($r);
    def second as Read init next($r);
    testing.assertEqual(hasNext($r), false);
    closeStrict($r);
    testing.assertEqual($first.id, "r1");
    testing.assertEqual($first.desc, "lane 1");
    testing.assertEqual($second.qual, "!!");
}

func testOpenStdinWrapsTheStdinSource() {
    # Opened and closed without reading: the test runner's stdin is not ours.
    # A stdin that cannot be reopened (a socket gives ENXIO) leaves nothing to
    # check here.
    def r as Reader;
    try {
        $r = openStdin();
    } catch (e) {
        return;
    }
    testing.assertEqual($r.src.origin, "<stdin>");
    testing.assertEqual($r.src.piped, false);
    testing.assertEqual(close($r), 0);
}
