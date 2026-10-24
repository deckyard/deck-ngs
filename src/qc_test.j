# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * White-box tests for the QC module.
 * Run with: jennifer test qc_test.j
 */
use testing;
use fs;
use strings;
use os;

func writeTemp(body as string, suffix as string) {
    def p as string init fs.makeTempFile("", "ngsqc", $suffix);
    fs.writeString($p, $body);
    return $p;
}

func fq(id as string, seq as string, qual as string) {
    return "@" + $id + "\n" + $seq + "\n+\n" + $qual + "\n";
}

func testCountsMatchHandComputedValues() {
    def body as string init fq("a", "ACGT", "IIII") + fq("b", "GGCC", "####");
    def p as string init writeTemp($body, ".fastq");
    def rep as Report init analyze($p, options());
    fs.remove($p);
    testing.assertEqual($rep.reads, 2);
    testing.assertEqual($rep.sampled, 2);
    testing.assertEqual($rep.bases, 8);
    testing.assertEqual($rep.gcBases, 6);
    testing.assertEqual($rep.minLength, 4);
    testing.assertEqual($rep.maxLength, 4);
    testing.assertEqual($rep.qualHist[40], 4);
    testing.assertEqual($rep.qualHist[2], 4);
    testing.assertEqual(meanQuality($rep), 21.0);
    testing.assertEqual(gcPercent($rep), 75.0);
    testing.assertEqual(q30Percent($rep), 50.0);
    testing.assertEqual(meanLength($rep), 4.0);
}

func testCountsAmbiguousBases() {
    def p as string init writeTemp(fq("a", "ACNNT", "IIIII"), ".fastq");
    def rep as Report init analyze($p, options());
    fs.remove($p);
    testing.assertEqual($rep.nBases, 2);
}

func testLengthHistogram() {
    def body as string init fq("a", "AC", "II") + fq("b", "ACGT", "IIII") + fq("c", "AC", "II");
    def p as string init writeTemp($body, ".fastq");
    def rep as Report init analyze($p, options());
    fs.remove($p);
    testing.assertEqual($rep.lengthHist[2], 2);
    testing.assertEqual($rep.lengthHist[4], 1);
    testing.assertEqual($rep.minLength, 2);
    testing.assertEqual($rep.maxLength, 4);
}

func testSubsamplingEveryNth() {
    def body as string init fq("a", "AA", "II") + fq("b", "CC", "II") + fq("c", "GG", "II") +
        fq("d", "TT", "II");
    def p as string init writeTemp($body, ".fastq");
    def opts as Options init options();
    $opts.everyNth = 2;
    def rep as Report init analyze($p, $opts);
    fs.remove($p);
    testing.assertEqual($rep.reads, 4);
    testing.assertEqual($rep.sampled, 2);
    testing.assertEqual($rep.bases, 4);
}

func testMaxReadsStopsEarly() {
    def body as string init fq("a", "AA", "II") + fq("b", "CC", "II") + fq("c", "GG", "II");
    def p as string init writeTemp($body, ".fastq");
    def opts as Options init options();
    $opts.maxReads = 2;
    def rep as Report init analyze($p, $opts);
    fs.remove($p);
    testing.assertEqual($rep.sampled, 2);
}

func testZeroValueOptionsAreUsable() {
    def p as string init writeTemp(fq("a", "ACGT", "IIII"), ".fastq");
    def bare as Options;
    def rep as Report init analyze($p, $bare);
    fs.remove($p);
    testing.assertEqual($rep.sampled, 1);
    testing.assertEqual(meanQuality($rep), 40.0);
}

func testPerBaseQuality() {
    def body as string init fq("a", "AC", "II") + fq("b", "GT", "#I");
    def p as string init writeTemp($body, ".fastq");
    def opts as Options init options();
    $opts.perBase = true;
    def rep as Report init analyze($p, $opts);
    fs.remove($p);
    def means as list of float init perBaseMean($rep);
    testing.assertEqual(len($means), 2);
    testing.assertEqual($means[0], 21.0);
    testing.assertEqual($means[1], 40.0);
}

func testPerBaseHandlesRaggedLengths() {
    def body as string init fq("a", "ACGT", "IIII") + fq("b", "AC", "II");
    def p as string init writeTemp($body, ".fastq");
    def opts as Options init options();
    $opts.perBase = true;
    def rep as Report init analyze($p, $opts);
    fs.remove($p);
    testing.assertEqual($rep.perBaseCount[0], 2);
    testing.assertEqual($rep.perBaseCount[3], 1);
}

func testDetectsAdapters() {
    def seq as string init "ACGT" + ADAPTER_TRUSEQ + "TT";
    def qual as string init strings.repeat("I", len($seq));
    def p as string init writeTemp(fq("a", $seq, $qual), ".fastq");
    def rep as Report init analyze($p, options());
    fs.remove($p);
    testing.assertEqual($rep.adapterHits[ADAPTER_TRUSEQ], 1);
    testing.assertEqual($rep.adapterHits[ADAPTER_NEXTERA], 0);
}

func testAdapterNotMatchedAcrossRecords() {
    def head as string init strings.substring(ADAPTER_TRUSEQ, 0, 6);
    def tail as string init strings.substring(ADAPTER_TRUSEQ, 6);
    def body as string init fq("a", $head, strings.repeat("I", len($head))) +
        fq("b", $tail, strings.repeat("I", len($tail)));
    def p as string init writeTemp($body, ".fastq");
    def rep as Report init analyze($p, options());
    fs.remove($p);
    testing.assertEqual($rep.adapterHits[ADAPTER_TRUSEQ], 0);
}

func testBatchBoundaryDoesNotChangeResults() {
    def body as string init "";
    for (def i in 0..25) {
        $body = $body + fq("r", "ACGTAC", "IIIII#");
    }
    def p as string init writeTemp($body, ".fastq");
    def one as Options init options();
    $one.batchSize = 1;
    def big as Options init options();
    $big.batchSize = 1000;
    def a as Report init analyze($p, $one);
    def b as Report init analyze($p, $big);
    fs.remove($p);
    testing.assertEqual($a.bases, $b.bases);
    testing.assertEqual($a.gcBases, $b.gcBases);
    testing.assertEqual(meanQuality($a), meanQuality($b));
}

func testGzippedInputGivesSameNumbers() {
    def body as string init fq("a", "ACGT", "IIII") + fq("b", "GGCC", "####");
    def p as string init writeTemp($body, ".fastq");
    def plain as Report init analyze($p, options());
    def gz as string init $p + ".gz";
    def zipped as os.Result init os.run(["gzip", "-kf", $p]);
    testing.assertEqual($zipped.exitCode, 0);
    def viaGzip as Report init analyze($gz, options());
    fs.remove($p);
    fs.remove($gz);
    testing.assertEqual($viaGzip.reads, $plain.reads);
    testing.assertEqual($viaGzip.bases, $plain.bases);
    testing.assertEqual($viaGzip.gcBases, $plain.gcBases);
    testing.assertEqual(meanQuality($viaGzip), meanQuality($plain));
}

func testEmptyFile() {
    def p as string init writeTemp("", ".fastq");
    def rep as Report init analyze($p, options());
    fs.remove($p);
    testing.assertEqual($rep.reads, 0);
    testing.assertEqual(meanQuality($rep), 0.0);
    testing.assertEqual(gcPercent($rep), 0.0);
    testing.assertEqual(meanLength($rep), 0.0);
}

func testWrongOffsetIsReported() {
    testing.assertThrows("analyzeWithBadOffset", "ngs.qc");
}

func analyzeWithBadOffset() {
    def p as string init writeTemp(fq("a", "ACGT", "!!!!"), ".fastq");
    def opts as Options init options();
    $opts.offset = fastq.PHRED64;
    def rep as Report init analyze($p, $opts);
    fs.remove($p);
    return $rep;
}

func testSummaryRenders() {
    def p as string init writeTemp(fq("a", "ACGT", "IIII"), ".fastq");
    def rep as Report init analyze($p, options());
    fs.remove($p);
    def text as string init summary($rep);
    testing.assertEqual(strings.contains($text, "reads seen"), true);
    testing.assertEqual(strings.contains($text, "GC content"), true);
}

func countDeckTemporaries() {
    def n as int init 0;
    for (def entry in fs.list(os.tempDir())) {
        if (strings.startsWith($entry, "ngsread") or strings.startsWith($entry, "ngswrite")) {
            $n = $n + 1;
        }
    }
    return $n;
}

func testFailedAnalyzeLeavesNoTemporaries() {
    # A malformed record part-way through a gzipped file makes the reader throw
    # while the decompressor and its FIFO are still live. Nothing may survive it.
    def body as string init fq("a", "ACGT", "IIII") + "@b\nACGT\nNOTAPLUS\nIIII\n";
    def p as string init writeTemp($body, ".fastq");
    def gz as string init $p + ".gz";
    def zipped as os.Result init os.run(["gzip", "-kf", $p]);
    testing.assertEqual($zipped.exitCode, 0);

    def before as int init countDeckTemporaries();
    def threw as bool init false;
    for (def attempt in 0..3) {
        try {
            def rep as Report init analyze($gz, options());
            testing.assertEqual($rep.reads, -1);
        } catch (e) {
            $threw = true;
        }
    }
    def after as int init countDeckTemporaries();
    fs.remove($p);
    fs.remove($gz);
    testing.assertEqual($threw, true);
    testing.assertEqual($after, $before);
}
