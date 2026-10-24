# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0
/**
 * White-box tests for the SAM reader.
 * Run with: jennifer test sam_test.j
 */
use testing;
use fs;

func temp(body as string, suffix as string) {
    def p as string init fs.makeTempFile("", "ngssam", $suffix);
    fs.writeString($p, $body);
    return $p;
}

func testParseAndFlags() {
    def head as string init "r1\t83\tchr1\t101\t60\t5M\t=\t200\t-150\t";
    def line as string init $head + "ACGTA\tIIIII\tNM:i:0\tAS:i:100";
    def rec as Record init parse($line);
    testing.assertEqual($rec.qname, "r1");
    testing.assertEqual($rec.pos, 101);
    testing.assertEqual(isReverse($rec), true);
    testing.assertEqual(isUnmapped($rec), false);
    testing.assertEqual(isPrimary($rec), true);
    testing.assertEqual(strandOf($rec), "-");
    testing.assertEqual(tagValue($rec, "NM", "?"), "0");
    testing.assertEqual(tagValue($rec, "XX", "?"), "?");
    testing.assertEqual(format($rec), $line);
}

func testCigar() {
    def ops as list of CigarOp init parseCigar("10M2I5D3S");
    testing.assertEqual(len($ops), 4);
    testing.assertEqual($ops[0].length, 10);
    testing.assertEqual($ops[1].op, "I");
    testing.assertEqual(referenceLength("10M2I5D3S"), 15);
    testing.assertEqual(queryLength("10M2I5D3S"), 15);
    testing.assertEqual(referenceLength("*"), 0);
}

func testToIntervalUsesCigarSpan() {
    def line as string init "r1\t0\tchr1\t101\t60\t10M5D10M\t*\t0\t0\t*\t*";
    def iv as intervals.Interval init toInterval(parse($line));
    testing.assertEqual($iv.start, 100);
    testing.assertEqual($iv.end, 125);
}

func testUnmappedIntervalFails() {
    testing.assertThrows("unmappedInterval", "ngs.sam");
}

func unmappedInterval() {
    return toInterval(parse("r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"));
}

func testHeaderAndReferenceSizes() {
    def head as string init "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:250\n";
    def body as string init $head + "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    def p as string init temp($body, ".sam");
    def r as Reader init open($p);
    def sizes as map of string to int init referenceSizes($r);
    def seen as list of string init [];
    def going as bool init true;
    while ($going) {
        def step as Next init next($r);
        match ($step) {
            when Alignment(hit) { $seen[] = $hit.value.qname; }
            when Done { $going = false; }
        }
    }
    closeStrict($r);
    fs.remove($p);
    testing.assertEqual(len($r.header), 3);
    testing.assertEqual($sizes["chr1"], 1000);
    testing.assertEqual($sizes["chr2"], 250);
    testing.assertEqual(len($seen), 1);
    testing.assertEqual($seen[0], "r1");
}

func flagged(flag as string) {
    return parse("r1\t" + $flag + "\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII");
}

func testFlagPredicates() {
    # 1040 = 0x400 duplicate + 0x10 reverse
    def dup as Record init flagged("1040");
    testing.assertEqual(hasFlag($dup, FLAG_DUPLICATE), true);
    testing.assertEqual(hasFlag($dup, FLAG_REVERSE), true);
    testing.assertEqual(hasFlag($dup, FLAG_PAIRED), false);
    testing.assertEqual(isDuplicate($dup), true);
    testing.assertEqual(isSecondary($dup), false);
    testing.assertEqual(isSupplementary($dup), false);
    testing.assertEqual(isPrimary($dup), true);

    def secondary as Record init flagged("256");
    testing.assertEqual(isSecondary($secondary), true);
    testing.assertEqual(isSupplementary($secondary), false);
    testing.assertEqual(isDuplicate($secondary), false);

    # 2064 = 0x800 supplementary + 0x10 reverse
    def supplementary as Record init flagged("2064");
    testing.assertEqual(isSupplementary($supplementary), true);
    testing.assertEqual(isSecondary($supplementary), false);
    testing.assertEqual(isPrimary($supplementary), false);

    # A combined mask asks whether any of its bits is set.
    testing.assertEqual(hasFlag($secondary, FLAG_SECONDARY | FLAG_SUPPLEMENTARY), true);
    testing.assertEqual(hasFlag(flagged("0"), FLAG_SECONDARY | FLAG_SUPPLEMENTARY), false);
}

func testOpenStdinLeavesTheHeaderEmpty() {
    # Opened and closed without reading: the test runner's stdin is not ours.
    # A stdin that cannot be reopened (a socket gives ENXIO) leaves nothing to
    # check here.
    def r as Reader;
    try {
        $r = openStdin();
    } catch (e) {
        return;
    }
    testing.assertEqual(len($r.header), 0);
    testing.assertEqual($r.src.origin, "<stdin>");
    testing.assertEqual(close($r), 0);
}
