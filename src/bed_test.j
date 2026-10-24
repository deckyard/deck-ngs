# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0
/**
 * White-box tests for the BED reader.
 * Run with: jennifer test bed_test.j
 */
use testing;
use fs;

func temp(body as string, suffix as string) {
    def p as string init fs.makeTempFile("", "ngsbed", $suffix);
    fs.writeString($p, $body);
    return $p;
}

func testParseAndFormat() {
    def iv as intervals.Interval init parse("chr1\t100\t200\tpeak1\t55\t-");
    testing.assertEqual($iv.start, 100);
    testing.assertEqual($iv.end, 200);
    testing.assertEqual($iv.name, "peak1");
    testing.assertEqual($iv.strand, "-");
    testing.assertEqual(format($iv), "chr1\t100\t200\tpeak1\t55\t-");
}

func testMinimalRecord() {
    def iv as intervals.Interval init parse("chr1\t0\t10");
    testing.assertEqual($iv.name, "");
    testing.assertEqual($iv.strand, ".");
    testing.assertEqual(format($iv), "chr1\t0\t10\t.\t0\t.");
}

func testSkipsDecorationLines() {
    def body as string init "# a comment\ntrack name=x\nchr1\t5\t10\nbrowser pos\nchr2\t1\t2\n";
    def p as string init temp($body, ".bed");
    def got as list of intervals.Interval init readAll($p);
    fs.remove($p);
    testing.assertEqual(len($got), 2);
    testing.assertEqual($got[0].chrom, "chr1");
    testing.assertEqual($got[1].chrom, "chr2");
}

func testTooFewColumnsFails() {
    testing.assertThrows("badBed", "ngs.bed");
}

func badBed() {
    return parse("chr1\t100");
}

func testWriteAllRoundTrip() {
    def p as string init fs.makeTempFile("", "ngsbed", ".bed");
    def xs as list of intervals.Interval init [intervals.named("chr1", 1, 5, "a", 10.0, "+")];
    writeAll($p, $xs);
    def back as list of intervals.Interval init readAll($p);
    fs.remove($p);
    testing.assertEqual(len($back), 1);
    testing.assertEqual($back[0].name, "a");
    testing.assertEqual($back[0].score, 10.0);
}

# New fixtures live inside the checkout, under the git-ignored .ci-fixtures/:
# in CI the interpreter runs in a container that mounts only the workspace.
func fixtureDir() {
    fs.mkdirAll(".ci-fixtures");
    return fs.makeTempDir(".ci-fixtures", "ngsbed");
}

func testWriteAppendsBed6Lines() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def p as string init $dir + "/out.bed";
    def k as source.Sink init source.createPlain($p);
    write($k, intervals.named("chr1", 100, 200, "peak1", 55.0, "+"));
    write($k, intervals.make("chr2", 0, 10));
    testing.assertEqual(source.closeSink($k), 0);
    testing.assertEqual(fs.readString($p), "chr1\t100\t200\tpeak1\t55\t+\nchr2\t0\t10\t.\t0\t.\n");
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
    testing.assertEqual(close($r), 0);
}
