# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * White-box tests for the trimming, filtering and subsampling module.
 * Run with: jennifer test trim_test.j
 */
use testing;
use fs;
use strings;
use math;

func writeTemp(body as string, suffix as string) {
    def p as string init fs.makeTempFile("", "ngstrim", $suffix);
    fs.writeString($p, $body);
    return $p;
}

func fq(id as string, seq as string, qual as string) {
    return "@" + $id + "\n" + $seq + "\n+\n" + $qual + "\n";
}

func testClipFixedEnds() {
    def rec as fastq.Read init fastq.parse("@r", "AACCGGTT", "+", "12345678");
    def got as fastq.Read init clipFixed($rec, 2, 2);
    testing.assertEqual($got.seq, "CCGG");
    testing.assertEqual($got.qual, "3456");
}

func testClipFixedOverlappingCutsGiveEmptyRead() {
    def rec as fastq.Read init fastq.parse("@r", "ACGT", "+", "IIII");
    def got as fastq.Read init clipFixed($rec, 3, 3);
    testing.assertEqual(len($got.seq), 0);
}

func testClipFixedNoOpWhenZero() {
    def rec as fastq.Read init fastq.parse("@r", "ACGT", "+", "IIII");
    testing.assertEqual(clipFixed($rec, 0, 0).seq, "ACGT");
}

func testAdapterClipFullMatch() {
    def seq as string init "ACGTACGT" + qc.ADAPTER_TRUSEQ;
    def rec as fastq.Read init fastq.parse("@r", $seq, "+", strings.repeat("I", len($seq)));
    def got as fastq.Read init clipAdapters($rec, qc.defaultAdapters(), 6);
    testing.assertEqual($got.seq, "ACGTACGT");
}

func testAdapterClipPartialAtReadEnd() {
    def prefix as string init strings.substring(qc.ADAPTER_TRUSEQ, 0, 8);
    def seq as string init "TTTTTTTT" + $prefix;
    def rec as fastq.Read init fastq.parse("@r", $seq, "+", strings.repeat("I", len($seq)));
    def got as fastq.Read init clipAdapters($rec, [qc.ADAPTER_TRUSEQ], 6);
    testing.assertEqual($got.seq, "TTTTTTTT");
}

func testAdapterClipRespectsMinOverlap() {
    def prefix as string init strings.substring(qc.ADAPTER_TRUSEQ, 0, 3);
    def seq as string init "TTTTTTTT" + $prefix;
    def rec as fastq.Read init fastq.parse("@r", $seq, "+", strings.repeat("I", len($seq)));
    testing.assertEqual(clipAdapters($rec, [qc.ADAPTER_TRUSEQ], 6).seq, $seq);
}

func testAdapterClipLeavesCleanRead() {
    def rec as fastq.Read init fastq.parse("@r", "ACGTACGTAC", "+", "IIIIIIIIII");
    testing.assertEqual(clipAdapters($rec, qc.defaultAdapters(), 6).seq, "ACGTACGTAC");
}

func testAdapterClipTakesEarliestOfSeveral() {
    def seq as string init "AAAA" + qc.ADAPTER_NEXTERA + "CCCC" + qc.ADAPTER_TRUSEQ;
    def rec as fastq.Read init fastq.parse("@r", $seq, "+", strings.repeat("I", len($seq)));
    testing.assertEqual(clipAdapters($rec, qc.defaultAdapters(), 6).seq, "AAAA");
}

func testQualityKeepLength() {
    testing.assertEqual(qualityKeepLength("IIIIIIII", 20, 2, 33), 8);
    testing.assertEqual(qualityKeepLength("IIII####", 20, 2, 33), 4);
    testing.assertEqual(qualityKeepLength("########", 20, 2, 33), 0);
    testing.assertEqual(qualityKeepLength("", 20, 2, 33), 0);
}

func testQualityTrimOnRecord() {
    def rec as fastq.Read init fastq.parse("@r", "ACGTACGT", "+", "IIII####");
    def got as fastq.Read init clipQuality($rec, 20, 2, 33);
    testing.assertEqual($got.seq, "ACGT");
    testing.assertEqual($got.qual, "IIII");
}

func testQualityTrimDisabledAtZeroThreshold() {
    def rec as fastq.Read init fastq.parse("@r", "ACGTACGT", "+", "########");
    testing.assertEqual(clipQuality($rec, 0, 2, 33).seq, "ACGTACGT");
}

func testTruncate() {
    def rec as fastq.Read init fastq.parse("@r", "ACGTACGT", "+", "IIIIIIII");
    testing.assertEqual(truncate($rec, 3).seq, "ACG");
    testing.assertEqual(truncate($rec, 0).seq, "ACGTACGT");
    testing.assertEqual(truncate($rec, 99).seq, "ACGTACGT");
}

func testKeepsFilters() {
    def opts as Options init options();
    $opts.minLength = 5;
    $opts.maxAmbiguous = 1;
    $opts.minMeanQuality = 30;
    testing.assertEqual(keeps(fastq.parse("@r", "ACGT", "+", "IIII"), $opts), false);
    testing.assertEqual(keeps(fastq.parse("@r", "ACGTAC", "+", "IIIIII"), $opts), true);
    testing.assertEqual(keeps(fastq.parse("@r", "ACNNAC", "+", "IIIIII"), $opts), false);
    testing.assertEqual(keeps(fastq.parse("@r", "ACGTAC", "+", "######"), $opts), false);
}

func testZeroValueOptionsKeepEverything() {
    def bare as Options;
    testing.assertEqual(keeps(fastq.parse("@r", "A", "+", "#"), $bare), true);
    testing.assertEqual(trimRead(fastq.parse("@r", "ACGT", "+", "IIII"), $bare).seq, "ACGT");
}

func testRunStreamsAndCounts() {
    def clean as string init "ACGTACGTACGTACGTACGTACGTACGT";
    def body as string init fq("keep", $clean, strings.repeat("I", len($clean))) +
        fq("short", "ACGT", "IIII") +
        fq(
        "adapter",
        $clean + qc.ADAPTER_TRUSEQ,
        strings.repeat("I", len($clean) + len(qc.ADAPTER_TRUSEQ)));
    def inPath as string init writeTemp($body, ".fastq");
    def outPath as string init fs.makeTempFile("", "ngsout", ".fastq");
    def st as Stats init run($inPath, $outPath, options());
    def written as string init fs.readString($outPath);
    fs.remove($inPath);
    fs.remove($outPath);
    testing.assertEqual($st.readsIn, 3);
    testing.assertEqual($st.readsOut, 2);
    testing.assertEqual($st.tooShort, 1);
    testing.assertEqual($st.adapterTrimmed, 1);
    testing.assertEqual(strings.contains($written, "@keep"), true);
    testing.assertEqual(strings.contains($written, "@short"), false);
    testing.assertEqual(strings.contains($written, qc.ADAPTER_TRUSEQ), false);
}

func testRunCountsBasesInAndOut() {
    def clean as string init "ACGTACGTACGTACGTACGTACGTACGT";
    def inPath as string init writeTemp(
        fq("a", $clean, strings.repeat("I", len($clean))),
        ".fastq");
    def outPath as string init fs.makeTempFile("", "ngsout", ".fastq");
    def opts as Options init options();
    $opts.trimFront = 3;
    def st as Stats init run($inPath, $outPath, $opts);
    fs.remove($inPath);
    fs.remove($outPath);
    testing.assertEqual($st.basesIn, len($clean));
    testing.assertEqual($st.basesOut, len($clean) - 3);
}

func testRunWritesGzipOutput() {
    def clean as string init "ACGTACGTACGTACGTACGTACGTACGT";
    def inPath as string init writeTemp(
        fq("a", $clean, strings.repeat("I", len($clean))),
        ".fastq");
    def outPath as string init fs.makeTempFile("", "ngsout", ".fastq.gz");
    fs.remove($outPath);
    def st as Stats init run($inPath, $outPath, options());
    def back as fastq.Reader init fastq.open($outPath);
    def n as int init 0;
    while (fastq.hasNext($back)) {
        def rec as fastq.Read init fastq.next($back);
        $n = $n + len($rec.id);
    }
    fastq.closeStrict($back);
    fs.remove($inPath);
    fs.remove($outPath);
    testing.assertEqual($st.readsOut, 1);
    testing.assertEqual($n, 1);
}

func testSubsampleEvery() {
    def body as string init "";
    for (def i in 0..10) {
        $body = $body + fq("r", "ACGT", "IIII");
    }
    def inPath as string init writeTemp($body, ".fastq");
    def outPath as string init fs.makeTempFile("", "ngssub", ".fastq");
    def st as Stats init subsampleEvery($inPath, $outPath, 3, 0);
    fs.remove($inPath);
    fs.remove($outPath);
    testing.assertEqual($st.readsIn, 10);
    testing.assertEqual($st.readsOut, 4);
}

func testSubsampleEveryRespectsMaxReads() {
    def body as string init "";
    for (def i in 0..10) {
        $body = $body + fq("r", "ACGT", "IIII");
    }
    def inPath as string init writeTemp($body, ".fastq");
    def outPath as string init fs.makeTempFile("", "ngssub", ".fastq");
    def st as Stats init subsampleEvery($inPath, $outPath, 1, 3);
    fs.remove($inPath);
    fs.remove($outPath);
    testing.assertEqual($st.readsOut, 3);
}

func testSubsampleFractionIsReproducibleWithSeed() {
    def body as string init "";
    for (def i in 0..50) {
        $body = $body + fq("r", "ACGT", "IIII");
    }
    def inPath as string init writeTemp($body, ".fastq");
    def outA as string init fs.makeTempFile("", "ngsfa", ".fastq");
    def outB as string init fs.makeTempFile("", "ngsfb", ".fastq");
    math.randSeed(42);
    def a as Stats init subsampleFraction($inPath, $outA, 0.5);
    math.randSeed(42);
    def b as Stats init subsampleFraction($inPath, $outB, 0.5);
    def sameBytes as bool init fs.readString($outA) == fs.readString($outB);
    fs.remove($inPath);
    fs.remove($outA);
    fs.remove($outB);
    testing.assertEqual($a.readsOut, $b.readsOut);
    testing.assertEqual($sameBytes, true);
}

func testSubsampleFractionRejectsBadRange() {
    testing.assertThrows("badFraction", "ngs.trim");
}

func badFraction() {
    return subsampleFraction("/dev/null", "/dev/null", 2.0);
}

func testRefusesToTrimAFileOntoItself() {
    # The output is opened for writing before the input is read, so allowing
    # this would truncate the input and report zero reads.
    def p as string init writeTemp(
        fq("a", "ACGTACGTACGTACGTACGTACGTACGT", strings.repeat("I", 28)),
        ".fastq");
    def sizeBefore as int init len(fs.readString($p));
    def refused as bool init false;
    try {
        def st as Stats init run($p, $p, options());
        testing.assertEqual($st.readsOut, -1);
    } catch (e) {
        $refused = true;
    }
    def sizeAfter as int init len(fs.readString($p));
    fs.remove($p);
    testing.assertEqual($refused, true);
    testing.assertEqual($sizeAfter, $sizeBefore);
}

func testRefusesToSubsampleAFileOntoItself() {
    def p as string init writeTemp(fq("a", "ACGT", "IIII"), ".fastq");
    def sizeBefore as int init len(fs.readString($p));
    def refused as bool init false;
    try {
        def st as Stats init subsampleEvery($p, $p, 1, 0);
        testing.assertEqual($st.readsOut, -1);
    } catch (e) {
        $refused = true;
    }
    def sizeAfter as int init len(fs.readString($p));
    fs.remove($p);
    testing.assertEqual($refused, true);
    testing.assertEqual($sizeAfter, $sizeBefore);
}
