# deck-ngs - `@jennifer/ngs`

[![test](https://github.com/jennifer-language/deck-ngs/actions/workflows/test.yml/badge.svg)](https://github.com/jennifer-language/deck-ngs/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-github%20pages-blue)](https://jennifer-language.github.io/deck-ngs/)
[![release](https://img.shields.io/github/v/release/jennifer-language/deck-ngs)](https://github.com/jennifer-language/deck-ngs/releases/latest)

Streaming file formats, read QC, interval arithmetic and pipeline orchestration
for next-generation sequencing data, in [Jennifer](https://jennifer-lang.dev/).

**[Read the documentation →](https://jennifer-language.github.io/deck-ngs/)**

## What it is

The glue and light-I/O layer of a sequencing pipeline, not the engine. Everything
streams, so a 50 GB input costs the same memory as a 50 MB one, and the heavy
compute stays in the native tools.

**In scope:** streaming format I/O (FASTQ, SAM, VCF, BED, GFF3/GTF), QC and
preprocessing, interval arithmetic, pipeline orchestration.

**Out of scope:** read alignment, assembly, variant calling and BAM/CRAM. Drive
those through the `pipeline` module and parse their text output here -
[why →](https://jennifer-language.github.io/deck-ngs/#out-of-scope).

## Install

```yaml
# deck.yaml
decks:
  "@jennifer/ngs": "^0.1.0"
```

```jennifer
import "@jennifer/ngs/" as ngs;            # the entry module
import "@jennifer/ngs/fastq.j" as fastq;   # the working surface
```

Needs **Jennifer 0.25.0 or newer**. Streaming gzip and pipeline orchestration
shell out, so they want the default `jennifer` binary rather than
`jennifer-tiny`, plus `gzip` and `mkfifo` on `PATH`; parsing, interval
arithmetic and formatting run on both. [Vendoring by hand, and the startup
capability check →](https://jennifer-language.github.io/deck-ngs/getting-started.html)

## A first look

```jennifer
# Streams, gzip and all - memory stays flat regardless of file size.
def report as qc.Report init qc.analyze("sample.fastq.gz", qc.options());
io.printf("%s", qc.summary($report));
```

## Modules

`source` · `fastq` · `bed` · `gff` · `sam` · `vcf` · `intervals` ·
`qc` · `trim` · `pipeline` · `core` · `ngs`

[What each one does →](https://jennifer-language.github.io/deck-ngs/#the-modules) - and every exported name, with its
signature, is in the [cheatsheet](https://jennifer-language.github.io/deck-ngs/cheatsheet.html).

## Documentation

The manual lives in `docs/` and is published to
[GitHub Pages](https://jennifer-language.github.io/deck-ngs/) by CI. Start with
[coordinate systems](https://jennifer-language.github.io/deck-ngs/coordinates.html) before writing interval code, and
[recipes](https://jennifer-language.github.io/deck-ngs/recipes.html) for worked end-to-end tasks.

Changing the deck - the gates, the docs build, CI and the release process - is
[Working on the deck](https://jennifer-language.github.io/deck-ngs/development.html).

## Licence

LGPL-3.0-only, matching the Jennifer interpreter. The full text is in
[`LICENSE`](LICENSE); every source file carries an SPDX header.
