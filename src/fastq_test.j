# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0
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

func testTrailingBlankLinesAreSkipped() {
    # LF and CRLF blank lines after the last record; hasNext used to report
    # them as a record and next then failed with an empty header.
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def p as string init $dir + "/blank.fastq";
    fs.writeString($p, "@a\nACGT\n+\nIIII\n\n\r\n\n");
    def r as Reader init open($p);
    def ids as list of string init [];
    while (hasNext($r)) {
        $ids[] = next($r).id;
    }
    closeStrict($r);
    testing.assertEqual($ids, ["a"]);
}

func testBlankLinesBeforeAHeaderAreSkipped() {
    # next without a hasNext probe skips them itself.
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def p as string init $dir + "/gap.fastq";
    fs.writeString($p, "@a\nAC\n+\nII\n\n@b\nGG\n+\n!!\n");
    def r as Reader init open($p);
    next($r);
    def second as Read init next($r);
    closeStrict($r);
    testing.assertEqual($second.id, "b");
}

func testNextPastTheEndRaisesNgsFastq() {
    # The header read sat outside the try, so this was a bare "runtime".
    testing.assertThrows("pastTheEnd", "ngs.fastq");
}

func pastTheEnd() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def p as string init $dir + "/one.fastq";
    fs.writeString($p, "@a\nAC\n+\nII\n");
    def r as Reader init open($p);
    next($r);
    return next($r).id;
}

func parseQuality(seq as string, qual as string) {
    return parse("@r", $seq, "+", $qual).id;
}

func testQualityOutsidePrintableAsciiIsRejected() {
    # 33..126 is the whole FASTQ alphabet. "é" is two bytes: it used to pass as
    # one rune and decode to two scores (195-33, 169-33).
    testing.assertThrowsWith("parseQuality", ["AC", "Ié"], "ngs.fastq");
    testing.assertThrowsWith("parseQuality", ["AC", "I "], "ngs.fastq");
    testing.assertThrowsWith("parseQuality", ["AC", "I\t"], "ngs.fastq");
    testing.assertThrowsWith("parseQuality", ["AC", "I\u007f"], "ngs.fastq");
    # Both ends of the range are accepted.
    testing.assertEqual(qualities(parse("@r", "AC", "+", "!~"), PHRED33), [0, 93]);
}

func testSequenceAndQualityLengthsCompareInBytes() {
    # "Aé" is 2 runes but 3 bytes against a 2-byte quality.
    testing.assertThrowsWith("parseQuality", ["Aé", "II"], "ngs.fastq");
}

func formatOf(id as string, desc as string, seq as string, qual as string) {
    return format(Read{id: $id, desc: $desc, seq: $seq, qual: $qual});
}

func testFormatRejectsFieldsThatBreakTheRecord() {
    # Each of these used to be written and then read back as something else.
    testing.assertThrowsWith("formatOf", ["a b", "", "A", "I"], "ngs.fastq");
    testing.assertThrowsWith("formatOf", ["x\nACGT\n+\nIIII\n@y", "", "A", "I"], "ngs.fastq");
    testing.assertThrowsWith("formatOf", ["a", "note\n@z", "A", "I"], "ngs.fastq");
    testing.assertThrowsWith("formatOf", ["a", "", "A\tC", "II"], "ngs.fastq");
    testing.assertThrowsWith("formatOf", ["a", "", "AC", "I"], "ngs.fastq");
    testing.assertThrowsWith("formatOf", ["a", "", "A", "\n"], "ngs.fastq");
    # A tab in the description is kept: samtools fastq -T writes tags there.
    testing.assertEqual(formatOf("a", "BC:Z:x\tRX:Z:y", "A", "I"), "@a BC:Z:x\tRX:Z:y\nA\n+\nI\n");
}
