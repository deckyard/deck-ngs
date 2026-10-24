# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0
/**
 * White-box tests for the GFF3 / GTF reader.
 * Run with: jennifer test gff_test.j
 */
use testing;
use fs;

func temp(body as string, suffix as string) {
    def p as string init fs.makeTempFile("", "ngsgff", $suffix);
    fs.writeString($p, $body);
    return $p;
}

func testGff3ParseAndCoordinateConversion() {
    def line as string init "chr1\tEnsembl\texon\t101\t200\t.\t+\t.\tID=e1;Name=Exon%3B1";
    def f as Feature init parse($line);
    testing.assertEqual($f.kind, "exon");
    testing.assertEqual($f.start, 101);
    testing.assertEqual($f.end, 200);
    testing.assertEqual(attribute($f, "ID", ""), "e1");
    testing.assertEqual(attribute($f, "Name", ""), "Exon;1");
    def iv as intervals.Interval init toInterval($f);
    testing.assertEqual($iv.start, 100);
    testing.assertEqual($iv.end, 200);
    testing.assertEqual(intervals.length($iv), 100);
}

func testGtfAttributeStyle() {
    def cols as string init "chr1\thavana\tCDS\t5\t9\t.\t-\t0\t";
    def line as string init $cols + "gene_id \"G1\"; transcript_id \"T1\";";
    def f as Feature init parse($line);
    testing.assertEqual(styleOf("gene_id \"G1\";"), STYLE_GTF);
    testing.assertEqual(attribute($f, "gene_id", ""), "G1");
    testing.assertEqual(identifier($f), "T1");
    testing.assertEqual($f.phase, 0);
}

func testRoundTripThroughInterval() {
    def iv as intervals.Interval init intervals.named("chr2", 10, 20, "f1", 0.0, "+");
    def f as Feature init fromInterval($iv, "ngs", "region");
    testing.assertEqual($f.start, 11);
    testing.assertEqual($f.end, 20);
    def back as intervals.Interval init toInterval($f);
    testing.assertEqual($back.start, 10);
    testing.assertEqual($back.end, 20);
}

func testAttributeEscaping() {
    def attrs as map of string to string init {"Note": "a;b=c,d"};
    def text as string init formatAttributes($attrs, STYLE_GFF3);
    testing.assertEqual($text, "Note=a%3Bb%3Dc%2Cd");
    def back as map of string to string init parseAttributes($text);
    testing.assertEqual($back["Note"], "a;b=c,d");
}

func testStopsAtFastaSection() {
    def gene as string init "chr1\t.\tgene\t1\t5\t.\t+\t.\tID=g1\n";
    def body as string init "##gff-version 3\n" + $gene + "##FASTA\n>chr1\nACGT\n";
    def p as string init temp($body, ".gff3");
    def got as list of Feature init readAll($p);
    fs.remove($p);
    testing.assertEqual(len($got), 1);
    testing.assertEqual(identifier($got[0]), "g1");
}

func testReadIntervalsFiltersByType() {
    def gene as string init "chr1\t.\tgene\t1\t5\t.\t+\t.\tID=g1\n";
    def exon as string init "chr1\t.\texon\t2\t4\t.\t+\t.\tID=e1\n";
    def body as string init $gene + $exon;
    def p as string init temp($body, ".gff3");
    def all as list of intervals.Interval init readIntervals($p, "");
    def exons as list of intervals.Interval init readIntervals($p, "exon");
    fs.remove($p);
    testing.assertEqual(len($all), 2);
    testing.assertEqual(len($exons), 1);
    testing.assertEqual($exons[0].start, 1);
}

func testUndefinedAttributeColumnRoundTrips() {
    # GFF3 writes "." for an undefined field. Parsing it must yield no
    # attributes, and re-formatting must write "." back - not ".=".
    def line as string init "chr1\t.\tgene\t1\t5\t.\t+\t.\t.";
    def f as Feature init parse($line);
    testing.assertEqual(len($f.attributes), 0);
    testing.assertEqual(identifier($f), "");
    testing.assertEqual(format($f, STYLE_GFF3), $line);
}

func testEmptyAttributesRenderAsDot() {
    def none as map of string to string init {};
    testing.assertEqual(formatAttributes($none, STYLE_GFF3), ".");
    testing.assertEqual(formatAttributes($none, STYLE_GTF), ".");
    testing.assertEqual(len(parseAttributes(".")), 0);
    testing.assertEqual(len(parseAttributes("")), 0);
}

func testGtfUnquotedNumericAttributes() {
    # GENCODE writes numeric GTF attributes unquoted: `exon_number 1; level 2;`
    def text as string init "gene_id \"G1\"; exon_number 1; level 2;";
    def m as map of string to string init parseAttributes($text);
    testing.assertEqual($m["gene_id"], "G1");
    testing.assertEqual($m["exon_number"], "1");
    testing.assertEqual($m["level"], "2");
    # ... and still detects GTF when the first pair itself is unquoted
    def bare as map of string to string init parseAttributes("gene_id G1; level 2;");
    testing.assertEqual($bare["gene_id"], "G1");
    testing.assertEqual($bare["level"], "2");
}

func testLeavesSpacesUnescaped() {
    # GFF3 allows unescaped spaces in attribute values and forbids encoding
    # anything outside the reserved set, so a space must survive untouched.
    def attrs as map of string to string init {"Note": "two words", "ID": "a;b"};
    def text as string init formatAttributes($attrs, STYLE_GFF3);
    testing.assertEqual(strings.contains($text, "two words"), true);
    testing.assertEqual(strings.contains($text, "%20"), false);
    testing.assertEqual(strings.contains($text, "%3B"), true);
    def back as map of string to string init parseAttributes($text);
    testing.assertEqual($back["Note"], "two words");
    testing.assertEqual($back["ID"], "a;b");
}

func testNamePrefersNameThenGeneNameThenIdentifier() {
    def both as Feature init parse("chr1\t.\tgene\t1\t9\t.\t+\t.\tID=g1;Name=BRCA2");
    testing.assertEqual(name($both), "BRCA2");
    def gtf as Feature init parse(
        "chr1\t.\texon\t1\t9\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"TP53\";");
    testing.assertEqual(name($gtf), "TP53");
    def idOnly as Feature init parse("chr1\t.\tgene\t1\t9\t.\t+\t.\tID=g7");
    testing.assertEqual(name($idOnly), "g7");
    def tx as string init "gene_id \"G1\"; transcript_id \"T1\";";
    def txOnly as Feature init parse("chr1\t.\texon\t1\t9\t.\t+\t.\t" + $tx);
    testing.assertEqual(name($txOnly), "T1");
    def bare as Feature init parse("chr1\t.\tregion\t1\t9\t.\t+\t.\t.");
    testing.assertEqual(name($bare), "");
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
