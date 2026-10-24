# deck-ngs - `@deckyard/ngs`

[![test](https://github.com/deckyard/deck-ngs/actions/workflows/test.yml/badge.svg)](https://github.com/deckyard/deck-ngs/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-github%20pages-blue)](https://deckyard.github.io/deck-ngs/)
[![release](https://img.shields.io/github/v/release/deckyard/deck-ngs)](https://github.com/deckyard/deck-ngs/releases/latest)

Streaming FASTQ, SAM, VCF, BED and GFF3/GTF I/O, read QC, trimming, interval
arithmetic and pipeline orchestration for next-generation sequencing data, in
[Jennifer](https://jennifer-lang.dev/). Memory use does not grow with file size.
Alignment, assembly, variant calling and BAM/CRAM stay with the native tools,
which the `pipeline` module drives.

## Install

```yaml
# deck.yaml
decks:
  "@deckyard/ngs": "^0.1.0"
```

`jvc` vendors the deck's `src/` into `vendor/deckyard/ngs/`. Release archives
are named `deckyard-ngs-<version>.tar.gz`;
[installing one by hand](https://deckyard.github.io/deck-ngs/getting-started.html)
is covered in the manual.

Needs Jennifer 0.25.0 or newer. Streaming gzip and `pipeline` need the `exec`
capability (the default `jennifer` binary, not `jennifer-tiny`) and `gzip` and
`mkfifo` on `PATH`; everything else runs on both binaries.

## Example

```jennifer
use io;
import "@deckyard/ngs/qc.j" as qc;

def report as qc.Report init qc.analyze("sample.fastq.gz", qc.options());
io.printf("%s", qc.summary($report));
```

## Documentation

The manual is at <https://deckyard.github.io/deck-ngs/>, built from `docs/`. The
[cheatsheet](https://deckyard.github.io/deck-ngs/cheatsheet.html) lists every
export.

## Licence

LGPL-3.0-only, see [LICENSE](LICENSE). Copyright (C) 2026 mplx &lt;mplx@jennifer-lang.dev&gt;.
