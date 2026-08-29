# @jennifer/ngs

Streaming file formats, read QC, interval arithmetic and pipeline orchestration
for next-generation sequencing data, written in [Jennifer](https://jennifer-lang.dev/).

## What this deck is

NGS breaks the ordinary deck model on two axes, and the design follows from
both.

**Scale.** FASTQ.gz files run to tens of gigabytes and BAM to hundreds. So
everything here streams: every reader is a pull cursor over one line at a time,
and nothing in this deck loads a sequencing file into memory. A 50 GB input
costs the same memory as a 50 MB one.

**Compute.** Alignment, assembly and variant calling are heavily SIMD-optimised
C. A tree-walking interpreter is orders of magnitude too slow to replace them,
and a library that pretended otherwise would be one nobody could use.

So the deck is deliberately **the glue and light-I/O layer, not the engine**.
It is the part of a sequencing pipeline that is all string handling, bookkeeping
and process management - the part where a real language with structs, enums,
error handling and a test runner beats three hundred lines of shell.

## In scope

| area | what you get |
| --- | --- |
| [Streaming format I/O](formats.md) | FASTQ, SAM, VCF, BED, GFF3/GTF - parse, filter and convert as a stream |
| [QC](qc.md) | per-base quality, length and GC distributions, adapter detection, one pass |
| [Preprocessing](preprocessing.md) | quality and adapter trimming, length filtering, subsampling |
| [Interval operations](intervals.md) | overlap, intersect, merge, subtract, coverage, complement |
| [Pipeline orchestration](pipelines.md) | drive the standard tools, fan out over samples, handle errors |

## Out of scope

Wrap the native tool and pipe to it; do not reimplement it here.

- **Read alignment** - BWA, Bowtie2, minimap2
- **De-novo assembly** - SPAdes
- **Variant calling** - GATK
- **BAM/CRAM handling** - additionally needs BGZF (blocked gzip for random
  access), which the `compress` library does not provide, and whose standalone
  value is limited without the downstream compute the deck omits

[Pipeline orchestration](pipelines.md) is how you drive all of them, and
[Recipes](recipes.md) shows what that looks like end to end.

## The modules

Every module is independent: importing `fastq.j` does not drag in the pipeline
layer. Import what you use.

| module | what it does |
| --- | --- |
| [`source`](streaming.md) | streaming line sources and sinks, gzip transparent both ways |
| [`fastq`](formats.md#fastq) | FASTQ records, streaming reader and writer, per-record accessors |
| [`bed`](formats.md#bed) | BED features, straight onto the interval engine |
| [`gff`](formats.md#gff3-and-gtf) | GFF3 and GTF features, with 1-based to 0-based conversion |
| [`sam`](formats.md#sam) | SAM alignments, FLAG and CIGAR handling |
| [`vcf`](formats.md#vcf) | VCF variants, INFO and per-sample FORMAT access |
| [`intervals`](intervals.md) | overlap, merge, intersect, subtract, coverage, complement, slop |
| [`qc`](qc.md) | one-pass FastQC-lite statistics |
| [`trim`](preprocessing.md) | fastp-lite trimming, filtering and subsampling |
| [`pipeline`](pipelines.md) | run, skip, check and fan out external tools |
| `core` | shared primitives the rest of the deck is built on |
| `ngs` | the entry module: version, and the host capability check |

Every exported name in all twelve is listed in the [cheatsheet](cheatsheet.md).

## A first look

```jennifer
use io;
import "@jennifer/ngs/fastq.j" as fastq;
import "@jennifer/ngs/qc.j" as qc;

# Streams, gzip and all - memory stays flat regardless of file size.
def report as qc.Report init qc.analyze("sample.fastq.gz", qc.options());
io.printf("%s", qc.summary($report));

def r as fastq.Reader init fastq.open("sample.fastq.gz");
while (fastq.hasNext($r)) {
    def rec as fastq.Read init fastq.next($r);
    io.printf("%s\t%d\t%f|prec=1\n", $rec.id, fastq.length($rec), fastq.gcPercent($rec));
}
fastq.closeStrict($r);
```

## Where to go next

- [Installing and importing](getting-started.md) - the vendor tree, and what
  your host needs to support
- [Streaming large files](streaming.md) - how gzip streaming works, and why it
  is built the way it is
- [Coordinate systems](coordinates.md) - the one page to read before writing
  interval code
- [Performance](performance.md) - what is fast, what is slow, and the measured
  numbers behind both

## Licence

LGPL-3.0-only, matching the Jennifer interpreter.
