# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

/**
 * White-box tests for the streaming source layer.
 * Run with: jennifer test source_test.j
 */
use testing;
use fs;
use os;
use convert;
use task;

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

# New fixtures live inside the checkout, under the git-ignored .ci-fixtures/:
# in CI the interpreter runs in a container that mounts only the workspace.
func fixtureDir() {
    fs.mkdirAll(".ci-fixtures");
    return fs.makeTempDir(".ci-fixtures", "ngssource");
}

func testOpenPlainStreamsLines() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def p as string init path.join($dir, "reads.fastq");
    fs.writeString($p, SAMPLE);
    def s as Source init openPlain($p);
    testing.assertEqual($s.piped, false);
    testing.assertEqual($s.helper, "");
    testing.assertEqual($s.origin, $p);
    def lines as list of string init [];
    while (hasNext($s)) {
        $lines[] = nextLine($s);
    }
    testing.assertEqual(close($s), 0);
    testing.assertEqual(len($lines), 8);
    testing.assertEqual($lines[0], "@r1");
    testing.assertEqual($lines[7], "))III");
}

func testOpenPlainRejectsAMissingFile() {
    testing.assertThrows("missingPlain", "ngs.source");
}

func missingPlain() {
    def s as Source init openPlain(".ci-fixtures/ngs-absent/none.fastq");
    return close($s);
}

func testOpenGzipRejectsAMissingFileBeforeSpawning() {
    testing.assertThrows("missingGzip", "ngs.source");
}

func missingGzip() {
    def s as Source init openGzip(".ci-fixtures/ngs-absent/none.fastq.gz");
    return close($s);
}

func testCreatePlainWritesVerbatim() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def p as string init path.join($dir, "out.txt");
    def k as Sink init createPlain($p);
    testing.assertEqual($k.piped, false);
    testing.assertEqual($k.helper, "");
    writeLine($k, "one");
    writeText($k, "two");
    testing.assertEqual(closeSink($k), 0);
    testing.assertEqual(fs.readString($p), "one\ntwo");
}

func testOpenStdinIsAPlainSource() {
    if (not fs.exists("/dev/stdin")) {
        testing.assertThrows("stdinSource", "ngs.source");
        return;
    }
    # Opened and closed without reading: the test runner's stdin is not ours.
    # A stdin that cannot be reopened (a socket gives ENXIO) leaves nothing to
    # check here.
    def s as Source;
    try {
        $s = openStdin();
    } catch (e) {
        return;
    }
    testing.assertEqual($s.origin, "<stdin>");
    testing.assertEqual($s.piped, false);
    testing.assertEqual($s.helper, "");
    testing.assertEqual(close($s), 0);
}

func stdinSource() {
    def s as Source init openStdin();
    return close($s);
}

# Runs `name` on a task and returns the kind it threw, "none" when it returned,
# or "timeout" when it was still blocked after 20 s. Used for calls whose
# regression was a hang: an open on a FIFO that never gets its other end.
func kindWithin(name as string, arg as string) {
    def t as task of string init spawn {
        try {
            if ($name == "createGzip") {
                closeSink(createGzip($arg));
            } elseif ($name == "open") {
                close(open($arg));
            } else {
                isGzip($arg);
            }
        } catch (e) {
            return $e.kind;
        }
        return "none";
    };
    try {
        return task.waitTimeout($t, 20000);
    } catch (e) {
        task.discard($t);
        return $e.kind;
    }
}

func testCreateGzipFailsFastOnAnUncreatableTarget() {
    # A regular file as the parent directory makes the target uncreatable
    # (ENOTDIR) for every user, root included. The old command opened the
    # target before the FIFO, died, and left our FIFO open waiting forever.
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def blocker as string init path.join($dir, "plain");
    fs.writeString($blocker, "x");
    testing.assertEqual(kindWithin("createGzip", path.join($blocker, "out.gz")), "ngs.source");
}

func testOpenRefusesAFifoInsteadOfLosingBytes() {
    # Sniffing used to read (and drop) the first two bytes of a FIFO, then
    # reopen it; with no writer the sniff itself blocked.
    # System temp, not .ci-fixtures: some shared mounts (vboxsf) refuse FIFOs.
    def dir as string init fs.makeTempDir("", "ngssrc");
    defer fs.removeAll($dir);
    def fifo as string init path.join($dir, "pipe");
    def made as os.Result init os.run(["mkfifo", $fifo]);
    testing.assertEqual($made.exitCode, 0);
    testing.assertEqual(kindWithin("open", $fifo), "ngs.source");
    testing.assertEqual(kindWithin("isGzip", $fifo), "ngs.source");
}

func testClosingAGzipSourceTwiceRaisesNgsSource() {
    testing.assertThrows("closeGzipTwice", "ngs.source");
}

func closeGzipTwice() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def p as string init path.join($dir, "reads.txt");
    fs.writeString($p, SAMPLE);
    def zipped as os.Result init os.run(["gzip", "-f", $p]);
    testing.assertEqual($zipped.exitCode, 0);
    def s as Source init open($p + ".gz");
    while (hasNext($s)) {
        nextLine($s);
    }
    closeStrict($s);
    return close($s);
}

func testClosingAGzipSinkTwiceRaisesNgsSource() {
    testing.assertThrows("closeGzipSinkTwice", "ngs.source");
}

func closeGzipSinkTwice() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def k as Sink init createGzip(path.join($dir, "out.txt.gz"));
    writeLine($k, "one");
    closeSink($k);
    return closeSink($k);
}
