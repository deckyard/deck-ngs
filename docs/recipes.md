# Recipes

Only the first recipe shows its imports; the others import each module they use
the same way. The runnable programs in `examples/` import from `../src/` instead.
Unless a recipe says otherwise, memory does not grow with the input.

## Convert FASTQ to FASTA

```jennifer
use io;
import "@deckyard/ngs/fastq.j" as fastq;
import "@deckyard/ngs/source.j" as source;

def r as fastq.Reader init fastq.open("in.fastq.gz");
def k as source.Sink init source.create("out.fasta.gz");
while (fastq.hasNext($r)) {
    source.writeText($k, fastq.toFasta(fastq.next($r)));
}
fastq.closeStrict($r);
source.closeSink($k);
```

## Filter reads by length and quality

For filters `trim.run` does not offer:

```jennifer
def r as fastq.Reader init fastq.open("in.fastq.gz");
def k as source.Sink init source.create("out.fastq.gz");
def kept as int init 0;

while (fastq.hasNext($r)) {
    def rec as fastq.Read init fastq.next($r);
    if (fastq.length($rec) < 100) { continue; }
    if (fastq.meanQuality($rec, fastq.PHRED33) < 25.0) { continue; }
    if (fastq.ambiguousCount($rec) > 2) { continue; }
    fastq.write($k, $rec);
    $kept = $kept + 1;
}
fastq.closeStrict($r);
source.closeSink($k);
io.printf("kept %d reads\n", $kept);
```

## QC a large run quickly

```jennifer
def opts as qc.Options init qc.options();
$opts.everyNth = 20;
$opts.maxReads = 200000;
def rep as qc.Report init qc.analyze("huge.fastq.gz", $opts);
io.printf("%s", qc.summary($rep));
```

## Trim adapters and low-quality tails

```jennifer
def opts as trim.Options init trim.options();
$opts.qualityThreshold = 20;
$opts.minLength = 36;
$opts.maxAmbiguous = 5;

def st as trim.Stats init trim.run("in.fastq.gz", "out.fastq.gz", $opts);
io.printf("%d/%d reads, %d/%d bases\n",
    $st.readsOut, $st.readsIn, $st.basesOut, $st.basesIn);
```

## Which peaks overlap exons?

```jennifer
def peaks as list of intervals.Interval init bed.readAll("peaks.bed");
def exons as list of intervals.Interval init gff.readIntervals("genes.gff3", "exon");

def sorted as list of intervals.Interval init intervals.sort($peaks);
def hits as list of int init intervals.countOverlaps($peaks, $exons);
def covered as list of int init intervals.coverage($peaks, $exons);

for (def i in 0..len($sorted)) {
    def iv as intervals.Interval init $sorted[$i];
    def fraction as float init $covered[$i] / intervals.length($iv);
    io.printf("%s\t%d\t%f|prec=3\n", intervals.toBed3($iv), $hits[$i], $fraction);
}
```

`countOverlaps` and `coverage` return results in the order of
`intervals.sort(peaks)`, not the input order, so the loop walks the sorted list.

## Which parts of my peaks are intergenic?

```jennifer
def genes as list of intervals.Interval init gff.readIntervals("genes.gff3", "gene");
def outside as list of intervals.Interval init intervals.subtract($peaks, $genes);
io.printf("%d intergenic pieces, %d bp\n",
    len($outside), intervals.totalLength($outside));
```

## Promoters from gene annotation

The 2 kb upstream of each gene, by strand:

```jennifer
def promoters as list of intervals.Interval init [];
for (def f in gff.readAll("genes.gff3")) {
    if ($f.kind != "gene") { continue; }
    def g as intervals.Interval init gff.toInterval($f);
    if ($g.strand == "-") {
        $promoters[] = intervals.named($g.chrom, $g.end, $g.end + 2000,
            gff.name($f), 0.0, "-");
    } else {
        $promoters[] = intervals.named($g.chrom, maxZero($g.start - 2000), $g.start,
            gff.name($f), 0.0, "+");
    }
}
bed.writeAll("promoters.bed", $promoters);

func maxZero(n as int) {
    if ($n < 0) { return 0; }
    return $n;
}
```

## Uncovered regions of the genome

```jennifer
def r as sam.Reader init sam.open("aln.sam");
def sizes as map of string to int init sam.referenceSizes($r);

def covered as list of intervals.Interval init [];
def going as bool init true;
while ($going) {
    def step as sam.Next init sam.next($r);
    match ($step) {
        when Alignment(hit) {
            if (not sam.isUnmapped($hit.value)) {
                $covered[] = sam.toInterval($hit.value);
            }
        }
        when Done { $going = false; }
    }
}
sam.closeStrict($r);

def gaps as list of intervals.Interval init intervals.complement($covered, $sizes);
bed.writeAll("uncovered.bed", $gaps);
```

`sam.toInterval` raises on an unmapped record, hence the filter. The list of
alignments has to fit in memory; for a real BAM, parse `samtools depth` output
instead.

## Read BAM through samtools

The deck reads SAM, not BAM. Piping `samtools view` into the program keeps the
SAM off the disk:

```sh
samtools view -h aln.bam | jennifer run myscript.j -
```

```jennifer
def r as sam.Reader init sam.openStdin();   # header is empty on stdin
```

Or run `samtools` from the program, through a file:

```jennifer
def j as pipeline.Job init pipeline.shell("extract",
    "samtools view -h -q 30 " + core.shellQuote($bam) + " > " + core.shellQuote($sam));
pipeline.check(pipeline.run($j));
def r as sam.Reader init sam.open($sam);
```

## High-confidence SNVs as BED

```jennifer
def r as vcf.Reader init vcf.open("calls.vcf.gz");
def out as source.Sink init source.create("snvs.bed");
def going as bool init true;
while ($going) {
    def step as vcf.Next init vcf.next($r);
    match ($step) {
        when Record(hit) { $going = emit($hit.value); }
        when Done { $going = false; }
    }
}
vcf.closeStrict($r);
source.closeSink($out);

func emit(v as vcf.Variant) {
    if (vcf.passed($v) and vcf.isSnv($v) and $v.qual >= 30.0) {
        bed.write($out, vcf.toInterval($v));
    }
    return true;
}
```

## Genotype counts across samples

```jennifer
def r as vcf.Reader init vcf.open("calls.vcf.gz");
def counts as map of string to int init {};
# ... walk records with vcf.next, then per variant:
for (def i in 0..len($r.samples)) {
    def call as string init vcf.genotype($v, $i);
    if (maps.has($counts, $call)) {
        $counts[$call] = $counts[$call] + 1;
    } else {
        $counts[$call] = 1;
    }
}
```

## A per-sample pipeline

The shape of `examples/pipeline.j`: trim in-process, then fan the external
step out across cores.

```jennifer
pipeline.requireTools(["bwa", "samtools"]);

def work as string init pipeline.workDir("run");
defer pipeline.cleanUp($work);

# Stage 1: trim each sample here.
def trimmed as list of string init [];
for (def sample in $samples) {
    def out as string init path.join($work, path.stem(path.base($sample)) + ".fq.gz");
    def st as trim.Stats init trim.run($sample, $out, trim.options());
    io.printf("%s: %d/%d reads\n", $sample, $st.readsOut, $st.readsIn);
    $trimmed[] = $out;
}

# Stage 2: hand the heavy work to the aligner, across cores.
def jobs as list of pipeline.Job init [];
for (def fq in $trimmed) {
    def bam as string init $fq + ".bam";
    def cmd as string init "bwa mem ref.fa " + core.shellQuote($fq)
        + " | samtools sort -o " + core.shellQuote($bam);
    def j as pipeline.Job init pipeline.shell(path.base($fq), $cmd);
    $jobs[] = pipeline.producing(pipeline.needing($j, [$fq]), [$bam]);
}

def results as list of pipeline.Result init pipeline.runParallel($jobs, 4);
io.printf("%s", pipeline.summary($results));
pipeline.checkAll($results);
```

Each job declares its inputs and outputs, so re-running after a crash skips the
jobs that finished. The trimming stage has no such check and runs again.
