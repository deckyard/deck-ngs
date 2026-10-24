# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0
/**
 * White-box tests for the VCF reader.
 * Run with: jennifer test vcf_test.j
 */
use testing;
use fs;

func temp(body as string, suffix as string) {
    def p as string init fs.makeTempFile("", "ngsvcf", $suffix);
    fs.writeString($p, $body);
    return $p;
}

func testParseAndInfo() {
    def head as string init "chr1\t500\trs1\tA\tG,T\t99.5\tPASS\tDP=30;AF=0.5;DB\t";
    def line as string init $head + "GT:DP\t0/1:20\t1/1:10";
    def v as Variant init parse($line);
    testing.assertEqual($v.chrom, "chr1");
    testing.assertEqual($v.pos, 500);
    testing.assertEqual($v.id, "rs1");
    testing.assertEqual(len($v.alt), 2);
    testing.assertEqual($v.qual, 99.5);
    testing.assertEqual(passed($v), true);
    testing.assertEqual(info($v, "DP", ""), "30");
    testing.assertEqual(hasInfo($v, "DB"), true);
    testing.assertEqual(info($v, "ZZ", "none"), "none");
    testing.assertEqual(genotype($v, 0), "0/1");
    testing.assertEqual(genotype($v, 1), "1/1");
    testing.assertEqual(sampleValue($v, 0, "DP", "?"), "20");
    testing.assertEqual(sampleValue($v, 5, "DP", "?"), "?");
}

func testVariantClassification() {
    def snv as Variant init parse("chr1\t1\t.\tA\tG\t.\t.\t.");
    def ins as Variant init parse("chr1\t1\t.\tA\tAGG\t.\t.\t.");
    def del as Variant init parse("chr1\t1\t.\tACGT\tA\t.\t.\t.");
    testing.assertEqual(isSnv($snv), true);
    testing.assertEqual(isIndel($snv), false);
    testing.assertEqual(isSnv($ins), false);
    testing.assertEqual(isIndel($ins), true);
    testing.assertEqual(isIndel($del), true);
}

func testSnvRequiresANucleotideReference() {
    # A padding or symbolic reference is not an SNV even when the alt is a base.
    testing.assertEqual(isSnv(parse("chr1\t1\t.\t-\tA\t.\t.\t.")), false);
    testing.assertEqual(isSnv(parse("chr1\t1\t.\t<DEL>\tA\t.\t.\t.")), false);
    testing.assertEqual(isSnv(parse("chr1\t1\t.\tA\tG\t.\t.\t.")), true);
    testing.assertEqual(isSnv(parse("chr1\t1\t.\tN\tG\t.\t.\t.")), true);
}

func testSymbolicAlleleIsNotIndel() {
    def sv as Variant init parse("chr1\t1\t.\tA\t<DEL>\t.\t.\t.");
    testing.assertEqual(isIndel($sv), false);
    testing.assertEqual(isSnv($sv), false);
}

func testIntervalSpansReferenceAllele() {
    def v as Variant init parse("chr1\t101\t.\tACGT\tA\t.\t.\t.");
    def iv as intervals.Interval init toInterval($v);
    testing.assertEqual($iv.start, 100);
    testing.assertEqual($iv.end, 104);
}

func testHeaderSamplesAndStream() {
    def meta as string init "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n";
    def cols as string init "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";
    def body as string init $meta + $cols + "chr1\t10\t.\tA\tG\t50\tPASS\tDP=9\tGT\t0/1\t1/1\n";
    def p as string init temp($body, ".vcf");
    def r as Reader init open($p);
    def n as int init 0;
    def going as bool init true;
    while ($going) {
        def step as Next init next($r);
        match ($step) {
            when Record(hit) { $n = $n + 1; }
            when Done { $going = false; }
        }
    }
    closeStrict($r);
    fs.remove($p);
    testing.assertEqual(len($r.meta), 2);
    testing.assertEqual(len($r.samples), 2);
    testing.assertEqual($r.samples[0], "S1");
    testing.assertEqual($r.samples[1], "S2");
    testing.assertEqual($n, 1);
}

func testFormatRoundTrip() {
    def line as string init "chr1\t500\trs1\tA\tG\t99.5\tPASS\tDP=30;DB\tGT\t0/1";
    testing.assertEqual(format(parse($line)), $line);
}

func testMinimalRecordRoundTrip() {
    def line as string init "chr1\t1\t.\tA\tG\t.\t.\t.";
    testing.assertEqual(format(parse($line)), $line);
}

func testParseInfo() {
    def m as map of string to string init parseInfo("DP=30;DB;AF=0.5;ANN=A|x=y");
    testing.assertEqual(len($m), 4);
    testing.assertEqual($m["DP"], "30");
    testing.assertEqual($m["DB"], "");
    testing.assertEqual($m["AF"], "0.5");
    # The value runs from the first `=` to the end of the entry.
    testing.assertEqual($m["ANN"], "A|x=y");
    testing.assertEqual(len(parseInfo(".")), 0);
    testing.assertEqual(len(parseInfo("")), 0);
    testing.assertEqual(len(parseInfo("DB;;H2")), 2);
}

func testFormatInfo() {
    testing.assertEqual(formatInfo({"DP": "30", "DB": "", "AF": "0.5"}), "DP=30;DB;AF=0.5");
    testing.assertEqual(formatInfo({}), ".");
    testing.assertEqual(formatInfo(parseInfo("NS=3;H2;AA=T")), "NS=3;H2;AA=T");
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
    testing.assertEqual(len($r.meta), 0);
    testing.assertEqual(len($r.samples), 0);
    testing.assertEqual(close($r), 0);
}
