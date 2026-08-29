# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * `@jennifer/ngs` - streaming I/O, QC, interval arithmetic and pipeline
 * orchestration for next-generation sequencing data.
 *
 * This entry module carries the deck's identity and a few whole-deck helpers.
 * The working surface lives in the submodules, imported directly:
 *
 * ```jennifer
 * import "@jennifer/ngs/fastq.j" as fastq;
 * import "@jennifer/ngs/qc.j" as qc;
 * import "@jennifer/ngs/intervals.j" as intervals;
 * import "@jennifer/ngs/pipeline.j" as pipeline;
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
 * **Scope.** The deck is the glue and light-I/O layer, not the heavy engine.
 * Alignment, assembly, variant calling and BAM/CRAM handling stay with the
 * native tools - drive them through `pipeline` and parse their text output here.
 *
 * @module ngs
 * @since 0.1.0
 */
use meta;
use strings;
import "./core.j" as core;

/** The deck's version. */
export def const VERSION as string init "0.1.0";

/** The deck's canonical package name. */
export def const NAME as string init "@jennifer/ngs";

/**
 * What this host can and cannot do with the deck.
 *
 * @field streamingGzip {bool} whether gzipped input can stream at constant
 *   memory; false falls back to `source.openGzipBuffered`
 * @field externalTools {bool} whether `pipeline` can run external programs
 * @field gzipOnPath {bool} whether a `gzip` executable was found on `PATH`
 */
export def struct Capabilities {
    streamingGzip as bool,
    externalTools as bool,
    gzipOnPath as bool
};

/**
 * What the current interpreter and host support.
 *
 * Streaming gzip and tool orchestration both need the `exec` capability, which
 * `jennifer-tiny` does not have. Check this once at startup rather than
 * discovering it partway through a long run.
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
 * A one-line description of the deck and what this host supports.
 *
 * @return {string} a human-readable status line
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
 * Raise unless the host can stream gzip and drive external tools.
 *
 * @throws {Error} on a host without the `exec` capability
 */
export func requireFullHost() {
    def caps as Capabilities init capabilities();
    if (not $caps.externalTools) {
        core.fail(
            "ngs",
            "this deck needs the `exec` capability for streaming gzip and tool orchestration");
    }
}
