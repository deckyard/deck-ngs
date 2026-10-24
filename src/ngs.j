# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * `@deckyard/ngs` - streaming I/O, QC, interval arithmetic and pipeline
 * orchestration for next-generation sequencing data.
 *
 * This entry module carries the deck's identity and the host capability check.
 * Import the submodules directly:
 *
 * ```jennifer
 * import "@deckyard/ngs/fastq.j" as fastq;
 * import "@deckyard/ngs/qc.j" as qc;
 * import "@deckyard/ngs/intervals.j" as intervals;
 * import "@deckyard/ngs/pipeline.j" as pipeline;
 * ```
 *
 * | module | what it does |
 * | --- | --- |
 * | `source` | streaming line sources and sinks, gzip transparent |
 * | `fastq` | FASTQ records, streaming reader and writer |
 * | `bed` | BED features, straight onto the interval engine |
 * | `gff` | GFF3 and GTF features, 1-based to 0-based conversion |
 * | `sam` | SAM alignments, FLAG and CIGAR handling |
 * | `vcf` | VCF variants, INFO and per-sample FORMAT access |
 * | `intervals` | overlap, merge, intersect, subtract, coverage |
 * | `qc` | one-pass FastQC-lite statistics |
 * | `trim` | fastp-lite trimming, filtering and subsampling |
 * | `pipeline` | run and fan out external tools |
 * | `core` | shared helpers the rest of the deck is built on |
 *
 * Alignment, assembly, variant calling and BAM/CRAM are out of scope: run the
 * native tools through `pipeline` and parse their text output here.
 *
 * @module ngs
 * @since 0.1.0
 */
use meta;
use strings;
import "./core.j" as core;

/** Deck version, kept by hand: `tools/inject-version.sh` does not stamp it. */
export def const VERSION as string init "0.1.0";

/** Package name. */
export def const NAME as string init "@deckyard/ngs";

/**
 * What this host can and cannot do with the deck.
 *
 * @field streamingGzip {bool} whether the host has `exec`, which streaming gzip
 *   needs; when false, use `source.openGzipBuffered`
 * @field externalTools {bool} whether `pipeline` can run external programs
 * @field gzipOnPath {bool} whether a `gzip` file was found on `PATH`
 */
export def struct Capabilities {
    streamingGzip as bool,
    externalTools as bool,
    gzipOnPath as bool
};

/**
 * What the current interpreter and host support.
 *
 * Streaming gzip and `pipeline` both need the `exec` capability, which
 * `jennifer-tiny` lacks. Check once at startup. `streamingGzip` does not
 * look at `gzipOnPath`; check both before streaming `.gz` files.
 *
 * @return {Capabilities} the host's capabilities
 */
export func capabilities() {
    def canExec as bool init meta.hasCapability("exec");
    return Capabilities{
        streamingGzip: $canExec,
        externalTools: $canExec,
        gzipOnPath: core.hasTool("gzip")
    };
}

/**
 * A one-line description of the deck and what this host supports, e.g.
 * `@deckyard/ngs 0.1.0 (streaming gzip, external tools)`.
 *
 * @return {string} the status line, no trailing newline
 */
export func describe() {
    def caps as Capabilities init capabilities();
    def notes as list of string init [];
    if ($caps.streamingGzip) {
        $notes[] = "streaming gzip";
    } else {
        $notes[] = "buffered gzip only";
    }
    if ($caps.externalTools) {
        $notes[] = "external tools";
    } else {
        $notes[] = "no external tools";
    }
    return NAME + " " + VERSION + " (" + strings.join($notes, ", ") + ")";
}

/**
 * Raise unless the host can stream gzip and run external tools. Does not
 * check for `gzip` on `PATH`.
 *
 * @throws {Error} kind "ngs.ngs" on a host without the `exec` capability
 */
export func requireFullHost() {
    def caps as Capabilities init capabilities();
    if (not $caps.externalTools) {
        core.fail(
            "ngs",
            "this deck needs the `exec` capability for streaming gzip and tool orchestration");
    }
}
