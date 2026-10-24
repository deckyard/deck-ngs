# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * White-box tests for the deck's shared primitives.
 * Run with: jennifer test core_test.j
 */
use testing;

func testCountIn() {
    testing.assertEqual(countIn("ACGTACGT", "A"), 2);
    testing.assertEqual(countIn("ACGTACGT", "ACGT"), 2);
    testing.assertEqual(countIn("AAAA", "AA"), 2);
    testing.assertEqual(countIn("ACGT", "X"), 0);
    testing.assertEqual(countIn("ACGT", ""), 0);
    testing.assertEqual(countIn("", "A"), 0);
}

func testGcAndNCounts() {
    testing.assertEqual(gcCount("ACGT"), 2);
    testing.assertEqual(gcCount("acgt"), 2);
    testing.assertEqual(gcCount("GGCC"), 4);
    testing.assertEqual(gcCount("AATT"), 0);
    testing.assertEqual(nCount("ACNNGT"), 2);
    testing.assertEqual(nCount("acnngt"), 2);
}

func testSplitTabAndFieldAccess() {
    def f as list of string init splitTab("a\tb\tc");
    testing.assertEqual(len($f), 3);
    testing.assertEqual(field($f, 1, "?"), "b");
    testing.assertEqual(field($f, 9, "?"), "?");
    testing.assertEqual(field($f, -1, "?"), "?");
}

func testTypedFieldAccess() {
    def f as list of string init splitTab("chr1\t100\t.\t2.5");
    testing.assertEqual(intField("test", $f, 1, -1), 100);
    testing.assertEqual(intField("test", $f, 2, -1), -1);
    testing.assertEqual(intField("test", $f, 9, -1), -1);
    testing.assertEqual(floatField("test", $f, 3, 0.0), 2.5);
    testing.assertEqual(floatField("test", $f, 2, 0.0), 0.0);
}

func testMalformedNumericFieldIsNamespaced() {
    testing.assertThrows("badNumericField", "ngs.test");
}

func badNumericField() {
    return intField("test", splitTab("chr1\tnotanint"), 1, 0);
}

func testParseHelpers() {
    testing.assertEqual(parseInt("test", "start", "42"), 42);
    testing.assertEqual(parseFloat("test", "qual", "2.5"), 2.5);
    testing.assertThrows("badParseInt", "ngs.test");
    testing.assertThrows("badParseFloat", "ngs.test");
}

func badParseInt() {
    return parseInt("test", "start", "zzz");
}

func badParseFloat() {
    return parseFloat("test", "qual", "zzz");
}

func testIsSkippable() {
    testing.assertEqual(isSkippable("", "#"), true);
    testing.assertEqual(isSkippable("   ", "#"), true);
    testing.assertEqual(isSkippable("# note", "#"), true);
    testing.assertEqual(isSkippable("@SQ\tSN:x", "@"), true);
    testing.assertEqual(isSkippable("chr1\t1\t2", "#"), false);
}

func testPadNum() {
    testing.assertEqual(padNum(5, 3), "005");
    testing.assertEqual(padNum(12345, 3), "12345");
    testing.assertEqual(padNum(0, 2), "00");
    testing.assertEqual(padNum(-7, 5), "-7");
}

func testPadNumHandlesWidthBeyondTheConstant() {
    # padNum is exported, so a width past its internal pad run must still pad
    # rather than fail inside its own slicing.
    testing.assertEqual(padNum(5, 25), "0000000000000000000000005");
    testing.assertEqual(len(padNum(7, 30)), 30);
}

func testPadNumOrdersNumerically() {
    testing.assertEqual(padNum(9, 4) < padNum(10, 4), true);
    testing.assertEqual(padNum(100, 4) < padNum(1000, 4), true);
}

func testShellQuote() {
    testing.assertEqual(shellQuote("plain"), "'plain'");
    testing.assertEqual(shellQuote("with space"), "'with space'");
    testing.assertEqual(shellQuote("a'b"), "'a'\\''b'");
    testing.assertEqual(shellQuote("; rm -rf /"), "'; rm -rf /'");
}

func testComplementBase() {
    testing.assertEqual(complementBase("A"), "T");
    testing.assertEqual(complementBase("g"), "c");
    testing.assertEqual(complementBase("N"), "N");
    testing.assertEqual(complementBase("R"), "Y");
    testing.assertEqual(complementBase("-"), "-");
}

func testReverseComplement() {
    testing.assertEqual(reverseComplement("ACGT"), "ACGT");
    testing.assertEqual(reverseComplement("AAAA"), "TTTT");
    testing.assertEqual(reverseComplement("ACGTN"), "NACGT");
    testing.assertEqual(reverseComplement("acgt"), "acgt");
    testing.assertEqual(reverseComplement(""), "");
}

func testReverseComplementIsAnInvolution() {
    def seq as string init "ACGTTGCANNRYACGT";
    testing.assertEqual(reverseComplement(reverseComplement($seq)), $seq);
}

func testHasTool() {
    testing.assertEqual(hasTool("sh"), true);
    testing.assertEqual(hasTool("definitelyNotARealToolXyz"), false);
}

func testFailUsesNamespacedKind() {
    testing.assertThrows("raiseSample", "ngs.sample");
}

func raiseSample() {
    return fail("sample", "boom");
}

func testReraiseWrapsAForeignError() {
    testing.assertThrows("reraiseRuntime", "ngs.bed");
    def message as string init "";
    try {
        reraiseRuntime();
    } catch (e) {
        $message = $e.message;
    }
    testing.assertEqual($message, "bad digit in record: chr1\tx\t5");
}

func reraiseRuntime() {
    return reraise("bed", "chr1\tx\t5", Error{kind: "runtime", message: "bad digit"});
}

func testReraiseKeepsADeckError() {
    testing.assertThrows("reraiseDeckError", "ngs.intervals");
}

func reraiseDeckError() {
    return reraise("bed", "chr1\t9\t5", Error{kind: "ngs.intervals", message: "inverted"});
}

# Fixtures live inside the checkout, under the git-ignored .ci-fixtures/: in CI
# the interpreter runs in a container that mounts only the workspace.
func fixtureDir() {
    fs.mkdirAll(".ci-fixtures");
    return fs.makeTempDir(".ci-fixtures", "ngscore");
}

func testSameFile() {
    def dir as string init fixtureDir();
    defer fs.removeAll($dir);
    def a as string init path.join($dir, "a.fastq");
    def b as string init path.join($dir, "b.fastq");
    fs.writeString($a, "@r\nA\n+\nI\n");
    fs.writeString($b, "@r\nA\n+\nI\n");
    testing.assertEqual(sameFile($a, $a), true);
    testing.assertEqual(sameFile($a, "./" + $a), true);
    testing.assertEqual(sameFile($a, $b), false);
    # A path that does not exist yet is never the same file, even as itself.
    def missing as string init path.join($dir, "out.fastq");
    testing.assertEqual(sameFile($a, $missing), false);
    testing.assertEqual(sameFile($missing, $missing), false);
}
