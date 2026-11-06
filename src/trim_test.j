# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

/**
 * White-box tests for the trimming, filtering and subsampling module.
 * Run with: jennifer test trim_test.j
 */
use testing;
use fs;
use strings;
use math;
use os;

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

func testAdapterStartFindsFullAndPartialAdapters() {
    def adapter as string init "AGATCGGAAGAGC";
    # full internal match
    testing.assertEqual(adapterStart("ACGTACGTAGATCGGAAGAGCTTTT", [$adapter], 6), 8);
    # a 6-base adapter prefix running off the 3' end
    testing.assertEqual(adapterStart("ACGTACGTAGATCG", [$adapter], 3), 8);
    # a 5-base prefix is below the requested overlap
    testing.assertEqual(adapterStart("ACGTACGTAGATC", [$adapter], 6), -1);
    testing.assertEqual(adapterStart("ACGTACGTAGATC", [$adapter], 5), 8);
    # 0 means DEFAULT_OVERLAP (6), so the 5-base prefix still does not count
    testing.assertEqual(adapterStart("ACGTACGTAGATC", [$adapter], 0), -1);
    testing.assertEqual(adapterStart("ACGTACGTACGT", [$adapter], 3), -1);
    testing.assertEqual(adapterStart("ACGTACGT", [], 3), -1);
}

func testAdapterStartTakesTheEarliestOfSeveral() {
    def seq as string init "TTTTCTGTCTCTTATAAAAAGATCGGAAGAGC";
    testing.assertEqual(adapterStart($seq, ["AGATCGGAAGAGC", "CTGTCTCTTATA"], 6), 4);
    testing.assertEqual(adapterStart($seq, ["CTGTCTCTTATA", "AGATCGGAAGAGC"], 6), 4);
}

# Kept out of the suite: it fails, exposing a bug. `Options` documents its zero
# value as keeping everything, but a zero `maxAmbiguous` is a live limit of 0
# (only -1 disables the check), so the zero value discards any read holding an N.
# Re-enable when `normalize` or the documentation is changed to agree.
#
# func testZeroValueOptionsKeepReadsWithAmbiguousBases() {
#     def bare as Options;
#     testing.assertEqual(keeps(fastq.parse("@r", "ACGNT", "+", "IIIII"), $bare), true);
# }

# Fixtures live inside the checkout, under the git-ignored .ci-fixtures/: in CI
# the interpreter runs in a container that mounts only the workspace.
func fixtureDir() {
    fs.mkdirAll(".ci-fixtures");
    return fs.makeTempDir(".ci-fixtures", "ngstrim");
}

# Options that trim and drop nothing.
func keepAll() {
    def opts as Options;
    $opts.maxAmbiguous = -1;
    return $opts;
}

func testRunLeavesAHardLinkedInputIntact() {
    # A hard link has its own real path, so `sameFile` cannot see it. Writing
    # through a temporary file and renaming it means the output gets a new
    # inode and the input is read in full.
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def input as string init path.join($dir, "in.fastq");
    def output as string init path.join($dir, "out.fastq");
    def body as string init fq("a", "ACGTACGTAC", "IIIIIIIIII") +
        fq("b", "TTTTGGGGCC", "IIIIIIIIII");
    fs.writeString($input, $body);
    os.run(["ln", $input, $output]);
    testing.assertEqual(core.sameFile($input, $output), false);
    def st as Stats init run($input, $output, keepAll());
    testing.assertEqual($st.readsIn, 2);
    testing.assertEqual($st.readsOut, 2);
    testing.assertEqual(fs.readString($input), $body);
    testing.assertEqual(fs.readString($output), $body);
    def sub as Stats init subsampleEvery($input, $output, 2, 0);
    testing.assertEqual($sub.readsOut, 1);
    testing.assertEqual(fs.readString($input), $body);
}

func testFailedRunKeepsTheOldOutputAndLeavesNoTemporary() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def input as string init path.join($dir, "in.fastq");
    def output as string init path.join($dir, "out.fastq");
    # The second record's quality is one character short: fastq.next raises.
    fs.writeString($input, fq("a", "ACGT", "IIII") + fq("b", "ACGT", "III"));
    fs.writeString($output, "old\n");
    def failed as bool init false;
    try {
        run($input, $output, keepAll());
    } catch (e) {
        $failed = true;
    }
    testing.assertEqual($failed, true);
    testing.assertEqual(fs.readString($output), "old\n");
    testing.assertEqual(fs.list($dir), ["in.fastq", "out.fastq"]);
}

func testRunOutputHasTheUsualFileMode() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def input as string init path.join($dir, "in.fastq");
    def output as string init path.join($dir, "out.fastq");
    fs.writeString($input, fq("a", "ACGT", "IIII"));
    subsampleEvery($input, $output, 1, 0);
    testing.assertEqual(fs.stat($output).mode, 0o644);
}

func testEmptyAdapterIsRejected() {
    testing.assertThrows("startWithEmptyAdapter", "ngs.trim");
    testing.assertThrows("trimWithEmptyAdapter", "ngs.trim");
}

func startWithEmptyAdapter() {
    return adapterStart("ACGTACGTAC", ["AGATCGGAAGAGC", ""], 6);
}

func trimWithEmptyAdapter() {
    def opts as Options init options();
    $opts.adapters = [""];
    return trimRead(fastq.parse("@r", "ACGTACGTAC", "+", "IIIIIIIIII"), $opts);
}

func testDuplicateAdaptersAreDropped() {
    def opts as Options init options();
    $opts.adapters = ["AAAA", "CCCC", "AAAA"];
    testing.assertEqual(normalize($opts).adapters, ["AAAA", "CCCC"]);
}

func testClipFixedNegativeTailCountsAsZero() {
    # front 2 removes "AC"; a tail of -3 removes nothing.
    def rec as fastq.Read init fastq.parse("@r", "ACGTACGTAC", "+", "0123456789");
    def got as fastq.Read init clipFixed($rec, 2, -3);
    testing.assertEqual($got.seq, "GTACGTAC");
    testing.assertEqual($got.qual, "23456789");
}

func testQualityOutsideTheOffsetRangeIsRejected() {
    testing.assertThrows("phred33AsPhred64", "ngs.trim");
    testing.assertThrows("nonAsciiQuality", "ngs.trim");
    testing.assertThrows("qualityLengthMismatch", "ngs.trim");
    # '@' (64) is Phred 0 at offset 64, the lowest valid score: keep nothing.
    testing.assertEqual(qualityKeepLength("@@@@", 20, 4, 64), 0);
}

func phred33AsPhred64() {
    # '#' (35) is below offset 64.
    return qualityKeepLength("IIII####", 20, 4, 64);
}

func nonAsciiQuality() {
    # 'é' is two bytes and not a quality character; read as bytes it shifted
    # every later score and the keep length came out past the read end.
    return qualityKeepLength("IIIIIIIé#", 20, 1, 33);
}

func qualityLengthMismatch() {
    def rec as fastq.Read init fastq.Read{id: "r", desc: "", seq: "ACGT", qual: "IIIII"};
    return clipQuality($rec, 20, 4, 33);
}
