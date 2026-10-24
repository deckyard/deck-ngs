# @deckyard/ngs

Streaming file formats, read QC, trimming, interval arithmetic and pipeline
orchestration for next-generation sequencing data, written in
[Jennifer](https://jennifer-lang.dev/).

This deck is part of [Deckyard](https://github.com/deckyard), community decks for [Jennifer](https://jennifer-lang.dev/).

## Scope

Every reader is a pull cursor over one line at a time, and nothing in the deck
loads a sequencing file into memory: a 50 GB input costs the same memory as a
50 MB one.

A tree-walking interpreter is orders of magnitude too slow to replace
SIMD-optimised C, so the deck is the glue and light-I/O layer of a pipeline:
string handling, bookkeeping and process management.

| area | what you get |
| --- | --- |
| [Format I/O](formats.md) | FASTQ, SAM, VCF, BED, GFF3/GTF: parse, filter and convert as a stream |
| [QC](qc.md) | quality, length and GC distributions, adapter counts, in one pass |
| [Preprocessing](preprocessing.md) | adapter and quality trimming, length and content filters, subsampling |
| [Interval operations](intervals.md) | overlap, intersect, merge, subtract, coverage, complement |
| [Pipeline orchestration](pipelines.md) | run the standard tools, skip finished steps, fan out over samples |

## Out of scope

Drive these through [`pipeline`](pipelines.md) and parse their text output:

- read alignment (BWA, Bowtie2, minimap2)
- de-novo assembly (SPAdes)
- variant calling (GATK)
- BAM and CRAM. These also need BGZF (blocked gzip with random access), which
  the `compress` library does not provide; pipe `samtools view -h` output into
  the [SAM reader](formats.md#sam) instead.

## The modules

Each module is imported on its own; importing `fastq.j` does not load the
pipeline layer.

| module | what it does |
| --- | --- |
| [`source`](streaming.md) | streaming line sources and sinks, gzip in both directions |
| [`fastq`](formats.md#fastq) | FASTQ records, reader, writer, per-record accessors |
| [`bed`](formats.md#bed) | BED features as `intervals.Interval` |
| [`gff`](formats.md#gff3-and-gtf) | GFF3 and GTF features, 1-based to 0-based conversion |
| [`sam`](formats.md#sam) | SAM alignments, FLAG and CIGAR handling |
| [`vcf`](formats.md#vcf) | VCF variants, INFO and per-sample FORMAT access |
| [`intervals`](intervals.md) | overlap, merge, intersect, subtract, coverage, complement, slop |
| [`qc`](qc.md) | one-pass FastQC-style statistics |
| [`trim`](preprocessing.md) | fastp-style trimming, filtering and subsampling |
| [`pipeline`](pipelines.md) | run, skip, check and fan out external tools |
| `core` | errors, field parsing, Go-speed counting, shell quoting |
| `ngs` | the entry module: name, version, host capability check |

Every export is in the [cheatsheet](cheatsheet.md).

## Example

```jennifer
use io;
import "@deckyard/ngs/fastq.j" as fastq;
import "@deckyard/ngs/qc.j" as qc;

def report as qc.Report init qc.analyze("sample.fastq.gz", qc.options());
io.printf("%s", qc.summary($report));

def r as fastq.Reader init fastq.open("sample.fastq.gz");
while (fastq.hasNext($r)) {
    def rec as fastq.Read init fastq.next($r);
    io.printf("%s\t%d\t%f|prec=1\n", $rec.id, fastq.length($rec), fastq.gcPercent($rec));
}
fastq.closeStrict($r);
```

Read [Coordinate systems](coordinates.md) before writing interval code, and
[Installing and importing](getting-started.md) for what the host needs.

## Licence

LGPL-3.0-only.
