# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * White-box tests for the streaming source layer.
 * Run with: jennifer test source_test.j
 */
use testing;
use fs;
use os;
use convert;

def const SAMPLE as string init "@r1\nACGTN\n+\nIIII#\n@r2\nGGCCA\n+\n))III\n";

func tempPath(suffix as string) {
    def p as string init fs.makeTempFile("", "ngstest", $suffix);
    fs.remove($p);
    return $p;
}

func testPlainRoundTrip() {
    def p as string init tempPath(".fastq");
    fs.writeString($p, SAMPLE);
    def s as Source init open($p);
    def n as int init 0;
    def first as string init "";
    while (hasNext($s)) {
        def line as string init nextLine($s);
        if ($n == 0) {
            $first = $line;
        }
        $n = $n + 1;
    }
    testing.assertEqual(close($s), 0);
    fs.remove($p);
    testing.assertEqual($n, 8);
    testing.assertEqual($first, "@r1");
}

func testGzipDetectionAndStreaming() {
    def p as string init tempPath(".fastq");
    fs.writeString($p, SAMPLE);
    def gz as string init $p + ".gz";
    def r as os.Result init os.run(["gzip", "-kf", $p]);
    testing.assertEqual($r.exitCode, 0);

    testing.assertEqual(isGzip($p), false);
    testing.assertEqual(isGzip($gz), true);

    def s as Source init open($gz);
    testing.assertEqual($s.piped, true);
    def lines as list of string init [];
    while (hasNext($s)) {
        $lines[] = nextLine($s);
    }
    closeStrict($s);
    testing.assertEqual(len($lines), 8);
    testing.assertEqual($lines[1], "ACGTN");
    testing.assertEqual($lines[7], "))III");
    fs.remove($p);
    fs.remove($gz);
}

func testGzipBufferedFallback() {
    def p as string init tempPath(".fastq");
    fs.writeString($p, SAMPLE);
    def gz as string init $p + ".gz";
    def r as os.Result init os.run(["gzip", "-kf", $p]);
    testing.assertEqual($r.exitCode, 0);

    def s as Source init openGzipBuffered($gz);
    def n as int init 0;
    while (hasNext($s)) {
        $n = $n + len(nextLine($s));
    }
    testing.assertEqual(close($s), 0);
    testing.assertEqual($n, 28);
    fs.remove($p);
    fs.remove($gz);
}

func testGzipSinkRoundTrip() {
    def out as string init tempPath(".txt.gz");
    def k as Sink init create($out);
    testing.assertEqual($k.piped, true);
    writeLine($k, "alpha");
    writeLine($k, "beta");
    testing.assertEqual(closeSink($k), 0);

    testing.assertEqual(isGzip($out), true);
    def s as Source init open($out);
    def got as list of string init [];
    while (hasNext($s)) {
        $got[] = nextLine($s);
    }
    closeStrict($s);
    testing.assertEqual(len($got), 2);
    testing.assertEqual($got[0], "alpha");
    testing.assertEqual($got[1], "beta");
    fs.remove($out);
}

func testPlainSinkRoundTrip() {
    def out as string init tempPath(".txt");
    def k as Sink init create($out);
    testing.assertEqual($k.piped, false);
    writeLine($k, "one");
    writeText($k, "two\n");
    testing.assertEqual(closeSink($k), 0);
    testing.assertEqual(fs.readString($out), "one\ntwo\n");
    fs.remove($out);
}

func testDoubleCloseReportsCleanlyAndStillCleansUp() {
    # A second close must raise the deck's own error rather than a bare fs one,
    # and a gzip source must still have released its FIFO and child either way.
    def p as string init tempPath(".fastq");
    fs.writeString($p, SAMPLE);
    def gz as string init $p + ".gz";
    def zipped as os.Result init os.run(["gzip", "-kf", $p]);
    testing.assertEqual($zipped.exitCode, 0);

    def s as Source init open($gz);
    def fifoDir as string init $s.helper;
    def drained as int init 0;
    while (hasNext($s)) {
        $drained = $drained + len(nextLine($s));
    }
    testing.assertEqual($drained, 28);
    closeStrict($s);
    testing.assertEqual(fs.exists($fifoDir), false);

    testing.assertThrows("closeTwice", "ngs.source");
    fs.remove($p);
    fs.remove($gz);
}

func closeTwice() {
    def p as string init tempPath(".txt");
    fs.writeString($p, "one\n");
    def s as Source init open($p);
    def second as int init close($s) + close($s);
    fs.remove($p);
    return $second;
}

func testMissingFileFails() {
    testing.assertThrows("missingOpen", "ngs.source");
}

func missingOpen() {
    def s as Source init open("/nonexistent/ngs/deck/nope.fastq");
    return close($s);
}

func testFifoCleanup() {
    def p as string init tempPath(".fastq");
    fs.writeString($p, SAMPLE);
    def gz as string init $p + ".gz";
    def r as os.Result init os.run(["gzip", "-kf", $p]);
    testing.assertEqual($r.exitCode, 0);
    def s as Source init open($gz);
    def fifo as string init $s.helper;
    testing.assertEqual(fs.exists($fifo), true);
    def drained as int init 0;
    while (hasNext($s)) {
        $drained = $drained + len(nextLine($s));
    }
    testing.assertEqual($drained, 28);
    closeStrict($s);
    testing.assertEqual(fs.exists($fifo), false);
    fs.remove($p);
    fs.remove($gz);
}
