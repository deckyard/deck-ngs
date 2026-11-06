# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0
/**
 * White-box tests for the SAM reader.
 * Run with: jennifer test sam_test.j
 */
use testing;
use fs;
use os;
use task;
use path;

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

func samLine(cigar as string, seq as string) {
    return parse("r1\t0\tchr1\t100\t60\t" + $cigar + "\t*\t0\t0\t" + $seq + "\t*");
}

func testToIntervalWithoutReferenceSpanCoversOneBase() {
    # POS 100 is 0-based 99. A "*" CIGAR used to take len(seq) (here 4) and a
    # clip-only CIGAR the clipped bases; htslib's bam_endpos gives pos + 1.
    def star as intervals.Interval init toInterval(samLine("*", "ACGT"));
    testing.assertEqual([$star.start, $star.end], [99, 100]);
    def clipped as intervals.Interval init toInterval(samLine("10S", "ACGTACGTAC"));
    testing.assertEqual([$clipped.start, $clipped.end], [99, 100]);
    def inserted as intervals.Interval init toInterval(samLine("2S3I", "ACGTA"));
    testing.assertEqual([$inserted.start, $inserted.end], [99, 100]);
}

func testCigarSpansAgreeWithParseCigar() {
    # 3S10M2D5M1I4N2=1X1H1P: reference M10 D2 M5 N4 =2 X1 = 24;
    # query S3 M10 M5 I1 =2 X1 = 22. Cross-checked against the parsed ops.
    def cigar as string init "3S10M2D5M1I4N2=1X1H1P";
    testing.assertEqual(referenceLength($cigar), 24);
    testing.assertEqual(queryLength($cigar), 22);
    def ref as int init 0;
    def query as int init 0;
    for (def op in parseCigar($cigar)) {
        if (strings.contains("MDN=X", $op.op)) {
            $ref = $ref + $op.length;
        }
        if (strings.contains("MIS=X", $op.op)) {
            $query = $query + $op.length;
        }
    }
    testing.assertEqual($ref, 24);
    testing.assertEqual($query, 22);
    testing.assertEqual(referenceLength("150M"), 150);
    testing.assertEqual(queryLength("*"), 0);
}

func testMalformedCigarSpansRaiseNgsSam() {
    testing.assertThrowsWith("refLength", ["M"], "ngs.sam");
    testing.assertThrowsWith("refLength", ["10M5"], "ngs.sam");
    testing.assertThrowsWith("refLength", ["9999999999999999999M"], "ngs.sam");
    # Ten 18-digit lengths, about 1e19, overflow the int64 total.
    testing.assertThrowsWith("refLength", [strings.repeat("999999999999999999M", 10)], "ngs.sam");
}

func refLength(cigar as string) {
    return referenceLength($cigar);
}

func formatWith(qname as string, tag as string) {
    def rec as Record init samLine("4M", "ACGT");
    $rec.qname = $qname;
    $rec.tags = [$tag];
    return format($rec);
}

func testFormatRejectsTabsAndLineBreaks() {
    # A tab shifted every later column; a newline split the record in two.
    testing.assertThrowsWith("formatWith", ["r\t1", "NM:i:0"], "ngs.sam");
    testing.assertThrowsWith("formatWith", ["r1", "NM:i:0\nr2"], "ngs.sam");
    testing.assertThrowsWith("formatWith", ["r1\r", "NM:i:0"], "ngs.sam");
    testing.assertEqual(
        formatWith("r1", "NM:i:0"),
        "r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0");
}

func testOpenRefusesAFifoInsteadOfBlocking() {
    # headerOf sniffed and read the FIFO, then open read it again from a
    # second handle; with no writer the first open blocked.
    # System temp, not .ci-fixtures: some shared mounts (vboxsf) refuse FIFOs.
    def dir as string init fs.makeTempDir("", "ngssam");
    defer fs.removeAll($dir);
    def fifo as string init path.join($dir, "pipe");
    testing.assertEqual(os.run(["mkfifo", $fifo]).exitCode, 0);
    def t as task of string init spawn {
        try {
            close(open($fifo));
        } catch (e) {
            return $e.kind;
        }
        return "none";
    };
    def kind as string init "timeout";
    try {
        $kind = task.waitTimeout($t, 20000);
    } catch (e) {
        task.discard($t);
    }
    testing.assertEqual($kind, "ngs.source");
}
