# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming read QC over FASTQ - the FastQC-lite half of the deck.
 *
 * One pass over the file produces read and base counts, the length
 * distribution, GC and N content, the full quality-score histogram, adapter
 * occurrence counts and (optionally) per-position mean quality.
 *
 * **How it stays fast.** A tree-walking interpreter costs roughly a microsecond
 * per step, so touching every base in Jennifer is what makes naive QC slow. This
 * module instead buffers a batch of records and counts characters with
 * `core.countIn`, which is a `strings.replace` inside Go - one native scan per
 * character class per batch rather than one interpreter step per base. That is
 * about ten times faster than the per-base loop and gives identical numbers.
 *
 * Per-position quality is the one statistic that cannot be expressed as a
 * character count, so `perBase` is opt-in and pairs naturally with `maxReads`.
 *
 * ```jennifer
 * def opts as qc.Options init qc.options();
 * $opts.maxReads = 200000;
 * def rep as qc.Report init qc.analyze("sample.fastq.gz", $opts);
 * io.printf("%s", qc.summary($rep));
 * ```
 *
 * @module qc
 * @since 0.1.0
 */
use io;
use strings;
use convert;
use maps;
import "./core.j" as core;
import "./fastq.j" as fastq;

/** Highest Phred score tracked; covers the whole printable ASCII range. */
export def const MAX_PHRED as int init 94;

/** Quality score at or above which a base is counted towards Q30. */
export def const Q30 as int init 30;

/** Records buffered before each batch of Go-speed counting. */
export def const DEFAULT_BATCH as int init 2000;

/**
 * Bases buffered before a batch is counted, whichever limit is reached first.
 *
 * `DEFAULT_BATCH` alone bounds a batch by record count, which is only a memory
 * bound when reads are short. Long-read data breaks that assumption badly: 2000
 * nanopore reads of 30 kb are 60 MB of sequence plus 60 MB of quality, and the
 * joins on top of that pushed one measured run to 657 MB resident. Bounding the
 * batch by bases as well keeps memory flat across read lengths, and costs
 * nothing for short reads, where the record limit is still the one that fires.
 */
export def const DEFAULT_BATCH_BASES as int init 4000000;

/** Illumina TruSeq / NEBNext universal adapter prefix. */
export def const ADAPTER_TRUSEQ as string init "AGATCGGAAGAGC";

/** Nextera / Illumina transposase adapter prefix. */
export def const ADAPTER_NEXTERA as string init "CTGTCTCTTATACACATCT";

/** Illumina small-RNA 3' adapter prefix. */
export def const ADAPTER_SMALL_RNA as string init "TGGAATTCTCGG";

/** Poly-G run, the two-colour-chemistry artefact seen on NovaSeq / NextSeq. */
export def const ARTEFACT_POLY_G as string init "GGGGGGGGGG";

/**
 * How to run an analysis.
 *
 * The zero value is valid and means "every read, no per-position quality,
 * Phred+33"; `options` returns the same thing explicitly.
 *
 * @field maxReads {int} stop after sampling this many reads; 0 means no limit
 * @field everyNth {int} sample one read in every N; 0 or 1 means every read
 * @field perBase {bool} also collect per-position quality (one step per base)
 * @field offset {int} quality encoding, `fastq.PHRED33` or `fastq.PHRED64`;
 *   0 means Phred+33
 * @field batchSize {int} records buffered per counting batch; 0 means the default
 * @field maxBatchBases {int} bases buffered per counting batch, whichever limit
 *   is reached first; 0 means `DEFAULT_BATCH_BASES`. Bounds memory on long-read
 *   data, where a record count does not
 * @field adapters {list of string} sequences to count; empty means the built-in set
 */
export def struct Options {
    maxReads as int,
    everyNth as int,
    perBase as bool,
    offset as int,
    batchSize as int,
    maxBatchBases as int,
    adapters as list of string
};

/**
 * Everything one pass measured.
 *
 * @field reads {int} records seen in the file
 * @field sampled {int} records actually measured, after subsampling
 * @field bases {int} bases in the sampled records
 * @field minLength {int} shortest sampled read, 0 when nothing was sampled
 * @field maxLength {int} longest sampled read
 * @field gcBases {int} G and C bases in the sampled records
 * @field nBases {int} ambiguous (N) bases in the sampled records
 * @field qualHist {list of int} count of bases per Phred score, indexed 0..93
 * @field lengthHist {map of int to int} read length to number of reads
 * @field adapterHits {map of string to int} adapter sequence to occurrences
 * @field perBaseSum {list of int} total quality per position, `perBase` only
 * @field perBaseCount {list of int} bases seen per position, `perBase` only
 */
export def struct Report {
    reads as int,
    sampled as int,
    bases as int,
    minLength as int,
    maxLength as int,
    gcBases as int,
    nBases as int,
    qualHist as list of int,
    lengthHist as map of int to int,
    adapterHits as map of string to int,
    perBaseSum as list of int,
    perBaseCount as list of int
};

/**
 * Default analysis options: every read, no per-position quality, Phred+33.
 *
 * @return {Options} a fresh, explicitly populated option set
 */
export func options() {
    def adapters as list of string init [];
    return Options{
        maxReads: 0,
        everyNth: 1,
        perBase: false,
        offset: fastq.PHRED33,
        batchSize: DEFAULT_BATCH,
        maxBatchBases: DEFAULT_BATCH_BASES,
        adapters: $adapters
    };
}

/**
 * The adapters and artefacts counted when `Options.adapters` is empty.
 *
 * @return {list of string} the built-in probe sequences
 */
export func defaultAdapters() {
    return [ADAPTER_TRUSEQ, ADAPTER_NEXTERA, ADAPTER_SMALL_RNA, ARTEFACT_POLY_G];
}

/**
 * Run a single streaming pass over a FASTQ file.
 *
 * Per-record work is kept to plain local counters and buffered strings; the
 * character scanning happens once per batch. Nothing proportional to the file
 * size is ever held in memory.
 *
 * @param path {string} `.fastq`, `.fq` or gzipped equivalent
 * @param opts {Options} what to measure; the zero value is valid
 * @return {Report} the collected statistics
 * @throws {Error} when the file is missing or a record is malformed
 */
export func analyze(path as string, opts as Options) {
    def cfg as Options init normalize($opts);
    def r as fastq.Reader init fastq.open($path);
    errdefer quietClose($r);
    def tally as Tally init emptyTally($cfg);
    def lengthHist as map of int to int init {};
    def sums as list of int init [];
    def counts as list of int init [];
    def seqs as list of string init [];
    def quals as list of string init [];
    def pending as int init 0;
    def reads as int init 0;
    def sampled as int init 0;
    def shortest as int init 0;
    def longest as int init 0;

    while (fastq.hasNext($r)) {
        if ($cfg.maxReads > 0 and $sampled >= $cfg.maxReads) {
            break;
        }
        def rec as fastq.Read init fastq.next($r);
        $reads = $reads + 1;
        if (not sampleThis($reads, $cfg)) {
            continue;
        }
        $sampled = $sampled + 1;
        def width as int init len($rec.seq);
        if ($sampled == 1 or $width < $shortest) {
            $shortest = $width;
        }
        if ($width > $longest) {
            $longest = $width;
        }
        if (maps.has($lengthHist, $width)) {
            $lengthHist[$width] = $lengthHist[$width] + 1;
        } else {
            $lengthHist[$width] = 1;
        }
        if ($cfg.perBase) {
            def raw as bytes init convert.bytesFromString($rec.qual, "utf-8");
            while (len($sums) < len($raw)) {
                $sums[] = 0;
                $counts[] = 0;
            }
            for (def i in 0..len($raw)) {
                $sums[$i] = $sums[$i] + ($raw[$i] -$cfg.offset);
                $counts[$i] = $counts[$i] + 1;
            }
        }
        $seqs[] = $rec.seq;
        $quals[] = $rec.qual;
        $pending = $pending + $width;
        if (len($seqs) >= $cfg.batchSize or $pending >= $cfg.maxBatchBases) {
            $tally = addBatch($tally, $seqs, $quals, $cfg);
            $seqs = [];
            $quals = [];
            $pending = 0;
        }
    }
    fastq.close($r);
    $tally = addBatch($tally, $seqs, $quals, $cfg);
    return Report{
        reads: $reads,
        sampled: $sampled,
        bases: $tally.bases,
        minLength: $shortest,
        maxLength: $longest,
        gcBases: $tally.gcBases,
        nBases: $tally.nBases,
        qualHist: $tally.qualHist,
        lengthHist: $lengthHist,
        adapterHits: $tally.adapterHits,
        perBaseSum: $sums,
        perBaseCount: $counts
    };
}

/**
 * Mean Phred quality across every sampled base.
 *
 * @param rep {Report} a completed report
 * @return {float} mean quality, 0.0 when nothing was sampled
 */
export func meanQuality(rep as Report) {
    if ($rep.bases == 0) {
        return 0.0;
    }
    def weighted as int init 0;
    for (def q in 0..MAX_PHRED) {
        $weighted = $weighted + $q * $rep.qualHist[$q];
    }
    return $weighted / $rep.bases;
}

/**
 * Percentage of sampled bases at Q30 or better.
 *
 * @param rep {Report} a completed report
 * @return {float} percentage in [0, 100], 0.0 when nothing was sampled
 */
export func q30Percent(rep as Report) {
    if ($rep.bases == 0) {
        return 0.0;
    }
    def good as int init 0;
    for (def q in Q30..MAX_PHRED) {
        $good = $good + $rep.qualHist[$q];
    }
    return 100.0 * $good / $rep.bases;
}

/**
 * G+C content of the sampled bases as a percentage.
 *
 * @param rep {Report} a completed report
 * @return {float} percentage in [0, 100], 0.0 when nothing was sampled
 */
export func gcPercent(rep as Report) {
    if ($rep.bases == 0) {
        return 0.0;
    }
    return 100.0 * $rep.gcBases / $rep.bases;
}

/**
 * Mean sampled read length.
 *
 * @param rep {Report} a completed report
 * @return {float} mean length, 0.0 when nothing was sampled
 */
export func meanLength(rep as Report) {
    if ($rep.sampled == 0) {
        return 0.0;
    }
    return $rep.bases / $rep.sampled;
}

/**
 * Mean quality at each read position.
 *
 * @param rep {Report} a report collected with `perBase` enabled
 * @return {list of float} mean quality per position, empty when `perBase` was off
 */
export func perBaseMean(rep as Report) {
    def out as list of float init [];
    for (def i in 0..len($rep.perBaseSum)) {
        if ($rep.perBaseCount[$i] == 0) {
            $out[] = 0.0;
        } else {
            $out[] = $rep.perBaseSum[$i] / $rep.perBaseCount[$i];
        }
    }
    return $out;
}

/**
 * Render a report as a human-readable text summary.
 *
 * @param rep {Report} a completed report
 * @return {string} a multi-line report, newline-terminated
 */
export func summary(rep as Report) {
    def out as string init "";
    $out = $out + io.sprintf("reads seen      %d\n", $rep.reads);
    $out = $out + io.sprintf("reads sampled   %d\n", $rep.sampled);
    $out = $out + io.sprintf("bases           %d\n", $rep.bases);
    $out = $out +
        io.sprintf(
        "read length     %d-%d (mean %f|prec=1)\n",
        $rep.minLength,
        $rep.maxLength,
        meanLength($rep));
    $out = $out + io.sprintf("GC content      %f|prec=2 %%\n", gcPercent($rep));
    $out = $out + io.sprintf("N bases         %d\n", $rep.nBases);
    $out = $out + io.sprintf("mean quality    %f|prec=2\n", meanQuality($rep));
    $out = $out + io.sprintf("bases >= Q%d     %f|prec=2 %%\n", Q30, q30Percent($rep));
    for (def probe in $rep.adapterHits) {
        $out = $out + io.sprintf("adapter %s %d\n", $probe, $rep.adapterHits[$probe]);
    }
    return $out;
}

/**
 * The batch-level accumulators, kept apart from `Report` so the hot per-record
 * path never copies a whole report.
 *
 * @field bases {int} bases counted so far
 * @field gcBases {int} G and C bases counted so far
 * @field nBases {int} N bases counted so far
 * @field qualHist {list of int} quality histogram so far
 * @field adapterHits {map of string to int} adapter occurrences so far
 */
def struct Tally {
    bases as int,
    gcBases as int,
    nBases as int,
    qualHist as list of int,
    adapterHits as map of string to int
};

/**
 * Fill in the zero value's implicit defaults.
 *
 * @param opts {Options} caller-supplied options
 * @return {Options} the same options with 0 fields resolved to their defaults
 */
func normalize(opts as Options) {
    def cfg as Options init $opts;
    if ($cfg.everyNth < 1) {
        $cfg.everyNth = 1;
    }
    if ($cfg.offset == 0) {
        $cfg.offset = fastq.PHRED33;
    }
    if ($cfg.batchSize < 1) {
        $cfg.batchSize = DEFAULT_BATCH;
    }
    if ($cfg.maxBatchBases < 1) {
        $cfg.maxBatchBases = DEFAULT_BATCH_BASES;
    }
    if (len($cfg.adapters) == 0) {
        $cfg.adapters = defaultAdapters();
    }
    return $cfg;
}

/**
 * A zeroed tally with its histogram and adapter table pre-sized.
 *
 * @param cfg {Options} normalized options
 * @return {Tally} an empty tally ready to accumulate into
 */
func emptyTally(cfg as Options) {
    def hist as list of int init [];
    for (def i in 0..MAX_PHRED) {
        $hist[] = 0;
    }
    def hits as map of string to int init {};
    for (def probe in $cfg.adapters) {
        $hits[$probe] = 0;
    }
    return Tally{bases: 0, gcBases: 0, nBases: 0, qualHist: $hist, adapterHits: $hits};
}

/**
 * Whether the Nth record read should be sampled.
 *
 * @param ordinal {int} 1-based position of the record in the file
 * @param cfg {Options} normalized options
 * @return {bool}
 */
func sampleThis(ordinal as int, cfg as Options) {
    return ($ordinal - 1) % $cfg.everyNth == 0;
}

/**
 * Count a buffered batch of records with native string scans.
 *
 * Sequences are joined with newlines so an adapter cannot be matched across a
 * record boundary; quality strings are joined bare, so the character histogram
 * counts exactly the sampled bases.
 *
 * @param t {Tally} the tally so far
 * @param seqs {list of string} buffered sequences
 * @param quals {list of string} buffered quality strings
 * @param cfg {Options} normalized options
 * @return {Tally} the updated tally
 */
func addBatch(t as Tally, seqs as list of string, quals as list of string, cfg as Options) {
    def out as Tally init $t;
    if (len($seqs) == 0) {
        return $out;
    }
    def joinedSeq as string init strings.join($seqs, "\n");
    def joinedQual as string init strings.join($quals, "");
    $out.bases = $out.bases + len($joinedQual);
    $out.gcBases = $out.gcBases + core.gcCount($joinedSeq);
    $out.nBases = $out.nBases + core.nCount($joinedSeq);
    for (def probe in $cfg.adapters) {
        $out.adapterHits[$probe] = $out.adapterHits[$probe] + core.countIn($joinedSeq, $probe);
    }
    $out.qualHist = tallyQuality($out.qualHist, $joinedQual, $cfg.offset);
    return $out;
}

/**
 * Add a batch of quality characters to the score histogram.
 *
 * Scans one Phred score at a time and stops as soon as every character in the
 * batch is accounted for, so a typical Illumina batch costs about forty native
 * passes rather than the full printable-ASCII range.
 *
 * @param hist {list of int} the histogram so far
 * @param chunk {string} concatenated quality strings
 * @param offset {int} quality encoding offset
 * @return {list of int} the updated histogram
 * @throws {Error} when a character decodes to a score outside the tracked range
 */
func tallyQuality(hist as list of int, chunk as string, offset as int) {
    def out as list of int init $hist;
    def total as int init len($chunk);
    def counted as int init 0;
    for (def q in 0..MAX_PHRED) {
        if ($counted >= $total) {
            break;
        }
        def hits as int init core.countIn($chunk, convert.fromCodepoint($q + $offset));
        if ($hits > 0) {
            $out[$q] = $out[$q] + $hits;
            $counted = $counted + $hits;
        }
    }
    if ($counted < $total) {
        core.fail("qc", "quality string holds characters outside Phred 0-93; wrong offset?");
    }
    return $out;
}

/**
 * Close a reader on a path that is already failing.
 *
 * Cleanup that runs while an error is propagating must not raise one of its own:
 * a deferred throw supersedes the pending error, which would replace the real
 * diagnosis with a complaint about an already-closed handle.
 *
 * @param r {fastq.Reader} the reader to close
 * @return {int} the close result, or -1 when closing itself failed
 */
func quietClose(r as fastq.Reader) {
    try {
        return fastq.close($r);
    } catch (e) {
        return -1;
    }
}
