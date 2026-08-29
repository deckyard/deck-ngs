# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Read trimming, filtering and subsampling over FASTQ - the fastp-lite half of
 * the deck's preprocessing.
 *
 * Every entry point streams: `run` reads one record, transforms it and writes it
 * out, so a 50 GB input needs no more memory than a 50 MB one. Output is
 * gzipped automatically when the target name ends in `.gz`.
 *
 * The trimming steps run in the same order fastp uses, each one exposed on its
 * own so a pipeline can apply just the parts it wants:
 *
 * 1. fixed clipping from either end (`clipFixed`)
 * 2. adapter removal from the 3' end (`clipAdapters`)
 * 3. sliding-window quality trimming (`clipQuality`)
 * 4. length capping (`truncate`)
 *
 * Filtering (`keeps`) is then applied to the trimmed record: minimum length,
 * maximum ambiguous bases, minimum mean quality.
 *
 * @module trim
 * @since 0.1.0
 */
use strings;
use convert;
use math;
import "./core.j" as core;
import "./source.j" as source;
import "./fastq.j" as fastq;
import "./qc.j" as qc;

/** Default sliding window width for quality trimming. */
export def const DEFAULT_WINDOW as int init 4;

/** Shortest partial adapter match accepted at a read's 3' end. */
export def const DEFAULT_OVERLAP as int init 6;

/**
 * How to trim and filter.
 *
 * The zero value is a no-op: it trims nothing and keeps everything, so you can
 * enable one behaviour at a time.
 *
 * @field trimFront {int} bases to remove from the 5' end
 * @field trimTail {int} bases to remove from the 3' end
 * @field adapters {list of string} adapters to clip; empty disables adapter
 *   clipping, `qc.defaultAdapters()` enables the built-in set
 * @field minOverlap {int} shortest partial adapter accepted at the 3' end;
 *   0 means `DEFAULT_OVERLAP`
 * @field qualityThreshold {int} sliding-window mean quality below which the read
 *   is cut; 0 disables quality trimming
 * @field windowSize {int} sliding window width; 0 means `DEFAULT_WINDOW`
 * @field maxLength {int} truncate reads longer than this; 0 means no cap
 * @field minLength {int} discard reads shorter than this after trimming
 * @field maxAmbiguous {int} discard reads with more than this many N bases;
 *   -1 disables the check
 * @field minMeanQuality {int} discard reads whose mean quality is below this;
 *   0 disables the check
 * @field offset {int} quality encoding; 0 means `fastq.PHRED33`
 */
export def struct Options {
    trimFront as int,
    trimTail as int,
    adapters as list of string,
    minOverlap as int,
    qualityThreshold as int,
    windowSize as int,
    maxLength as int,
    minLength as int,
    maxAmbiguous as int,
    minMeanQuality as int,
    offset as int
};

/**
 * What a run did.
 *
 * @field readsIn {int} records read
 * @field readsOut {int} records written
 * @field basesIn {int} bases read
 * @field basesOut {int} bases written
 * @field adapterTrimmed {int} records shortened by adapter clipping
 * @field qualityTrimmed {int} records shortened by quality trimming
 * @field tooShort {int} records discarded for being shorter than `minLength`
 * @field tooAmbiguous {int} records discarded for exceeding `maxAmbiguous`
 * @field lowQuality {int} records discarded for failing `minMeanQuality`
 */
export def struct Stats {
    readsIn as int,
    readsOut as int,
    basesIn as int,
    basesOut as int,
    adapterTrimmed as int,
    qualityTrimmed as int,
    tooShort as int,
    tooAmbiguous as int,
    lowQuality as int
};

/**
 * Default options: adapter clipping on, sliding-window quality trimming at Q20,
 * and reads shorter than 25 bases discarded.
 *
 * @return {Options} a populated option set to adjust and pass to `run`
 */
export func options() {
    return Options{
        trimFront: 0,
        trimTail: 0,
        adapters: qc.defaultAdapters(),
        minOverlap: DEFAULT_OVERLAP,
        qualityThreshold: 20,
        windowSize: DEFAULT_WINDOW,
        maxLength: 0,
        minLength: 25,
        maxAmbiguous: -1,
        minMeanQuality: 0,
        offset: fastq.PHRED33
    };
}

/**
 * Stream a FASTQ file through trimming and filtering into a new file.
 *
 * @param inPath {string} input FASTQ, optionally gzipped
 * @param outPath {string} output FASTQ; gzipped when the name ends in `.gz`
 * @param opts {Options} what to trim and filter; the zero value is a no-op copy
 * @return {Stats} what the run did
 * @throws {Error} when either file cannot be opened, when the compressor fails,
 *   or when `inPath` and `outPath` are the same file
 */
export func run(inPath as string, outPath as string, opts as Options) {
    if (core.sameFile($inPath, $outPath)) {
        core.fail(
            "trim",
            "input and output are the same file (" + $outPath +
                "); the output is opened for writing before the input is read, " +
                "which would truncate it");
    }
    def cfg as Options init normalize($opts);
    def r as fastq.Reader init fastq.open($inPath);
    errdefer quietClose($r);
    def k as source.Sink init source.create($outPath);
    errdefer quietCloseSink($k);
    def st as Stats init Stats{
        readsIn: 0,
        readsOut: 0,
        basesIn: 0,
        basesOut: 0,
        adapterTrimmed: 0,
        qualityTrimmed: 0,
        tooShort: 0,
        tooAmbiguous: 0,
        lowQuality: 0
    };

    while (fastq.hasNext($r)) {
        def rec as fastq.Read init fastq.next($r);
        $st.readsIn = $st.readsIn + 1;
        $st.basesIn = $st.basesIn + len($rec.seq);

        def work as fastq.Read init clipFixed($rec, $cfg.trimFront, $cfg.trimTail);
        def afterFixed as int init len($work.seq);
        $work = clipAdapters($work, $cfg.adapters, $cfg.minOverlap);
        if (len($work.seq) < $afterFixed) {
            $st.adapterTrimmed = $st.adapterTrimmed + 1;
        }
        def afterAdapter as int init len($work.seq);
        $work = clipQuality($work, $cfg.qualityThreshold, $cfg.windowSize, $cfg.offset);
        if (len($work.seq) < $afterAdapter) {
            $st.qualityTrimmed = $st.qualityTrimmed + 1;
        }
        $work = truncate($work, $cfg.maxLength);

        # `judge` already applies every filter to record why a read was dropped,
        # so asking `keeps` the same question would re-run them - including the
        # ambiguous-base count and the per-base mean quality.
        def kept as int init $st.readsOut;
        $st = judge($st, $work, $cfg);
        if ($st.readsOut > $kept) {
            fastq.write($k, $work);
        }
    }
    fastq.close($r);
    def code as int init source.closeSink($k);
    if ($code != 0) {
        core.fail(
            "trim",
            "compressor for " + $outPath + " exited with code " + convert.toString($code));
    }
    return $st;
}

/**
 * Apply every configured trimming step to one record.
 *
 * @param rec {Read} the record to trim
 * @param opts {Options} what to trim; the zero value returns the record unchanged
 * @return {Read} the trimmed record, possibly empty
 */
export func trimRead(rec as fastq.Read, opts as Options) {
    def cfg as Options init normalize($opts);
    def work as fastq.Read init clipFixed($rec, $cfg.trimFront, $cfg.trimTail);
    $work = clipAdapters($work, $cfg.adapters, $cfg.minOverlap);
    $work = clipQuality($work, $cfg.qualityThreshold, $cfg.windowSize, $cfg.offset);
    return truncate($work, $cfg.maxLength);
}

/**
 * Remove a fixed number of bases from each end.
 *
 * @param rec {Read} the record
 * @param front {int} bases to remove from the 5' end
 * @param tail {int} bases to remove from the 3' end
 * @return {Read} the clipped record, empty when the cuts overlap
 */
export func clipFixed(rec as fastq.Read, front as int, tail as int) {
    if ($front <= 0 and $tail <= 0) {
        return $rec;
    }
    def start as int init $front;
    if ($start < 0) {
        $start = 0;
    }
    def stop as int init len($rec.seq) - $tail;
    if ($stop < $start) {
        $stop = $start;
    }
    if ($start > len($rec.seq)) {
        $start = len($rec.seq);
        $stop = $start;
    }
    return fastq.subseq($rec, $start, $stop);
}

/**
 * Clip the read at the first adapter occurrence, including a partial adapter
 * running off the 3' end.
 *
 * A full internal match is found with a native substring scan; a partial match
 * is found by testing adapter prefixes against the read's suffix, longest first,
 * down to `minOverlap` bases.
 *
 * @param rec {Read} the record
 * @param adapters {list of string} adapter sequences; empty disables clipping
 * @param minOverlap {int} shortest partial adapter accepted at the 3' end
 * @return {Read} the record clipped at the earliest adapter start
 */
export func clipAdapters(rec as fastq.Read, adapters as list of string, minOverlap as int) {
    if (len($adapters) == 0 or len($rec.seq) == 0) {
        return $rec;
    }
    def cut as int init adapterStart($rec.seq, $adapters, $minOverlap);
    if ($cut < 0) {
        return $rec;
    }
    return fastq.subseq($rec, 0, $cut);
}

/**
 * The earliest position at which an adapter starts in a sequence.
 *
 * @param seq {string} the sequence to scan
 * @param adapters {list of string} adapter sequences
 * @param minOverlap {int} shortest partial adapter accepted at the 3' end
 * @return {int} zero-based cut position, or -1 when no adapter is present
 */
export func adapterStart(seq as string, adapters as list of string, minOverlap as int) {
    def overlap as int init $minOverlap;
    if ($overlap < 1) {
        $overlap = DEFAULT_OVERLAP;
    }
    def best as int init -1;
    for (def probe in $adapters) {
        def at as int init strings.indexOf($seq, $probe);
        if ($at < 0) {
            $at = partialAt($seq, $probe, $overlap);
        }
        if ($at >= 0 and ($best < 0 or $at < $best)) {
            $best = $at;
        }
    }
    return $best;
}

/**
 * Trim the 3' end where a sliding window's mean quality falls below a threshold.
 *
 * This is the Trimmomatic / fastp sliding-window rule: scan windows from the 5'
 * end and cut at the first window whose mean quality is too low.
 *
 * @param rec {Read} the record
 * @param threshold {int} mean Phred quality a window must reach; 0 disables
 * @param window {int} window width; 0 means `DEFAULT_WINDOW`
 * @param offset {int} quality encoding offset
 * @return {Read} the trimmed record
 */
export func clipQuality(rec as fastq.Read, threshold as int, window as int, offset as int) {
    if ($threshold <= 0 or len($rec.qual) == 0) {
        return $rec;
    }
    def keep as int init qualityKeepLength($rec.qual, $threshold, $window, $offset);
    if ($keep >= len($rec.seq)) {
        return $rec;
    }
    return fastq.subseq($rec, 0, $keep);
}

/**
 * How many bases survive sliding-window quality trimming.
 *
 * @param qual {string} the raw quality string
 * @param threshold {int} mean Phred quality a window must reach
 * @param window {int} window width; 0 means `DEFAULT_WINDOW`
 * @param offset {int} quality encoding offset
 * @return {int} number of bases to keep from the 5' end
 */
export func qualityKeepLength(qual as string, threshold as int, window as int, offset as int) {
    def width as int init $window;
    if ($width < 1) {
        $width = DEFAULT_WINDOW;
    }
    def raw as bytes init convert.bytesFromString($qual, "utf-8");
    def n as int init len($raw);
    if ($n == 0) {
        return 0;
    }
    if ($width > $n) {
        $width = $n;
    }
    def running as int init 0;
    for (def i in 0..$width) {
        $running = $running + ($raw[$i] -$offset);
    }
    def limit as int init $threshold * $width;
    if ($running < $limit) {
        return 0;
    }
    for (def start in 1..($n - $width + 1)) {
        $running = $running - ($raw[$start - 1] -$offset) + ($raw[$start + $width - 1] -$offset);
        if ($running < $limit) {
            return $start;
        }
    }
    return $n;
}

/**
 * Cap a record's length, keeping the 5' end.
 *
 * @param rec {Read} the record
 * @param maxLength {int} maximum length; 0 or less means no cap
 * @return {Read} the capped record
 */
export func truncate(rec as fastq.Read, maxLength as int) {
    if ($maxLength <= 0 or len($rec.seq) <= $maxLength) {
        return $rec;
    }
    return fastq.subseq($rec, 0, $maxLength);
}

/**
 * Whether a trimmed record passes the configured filters.
 *
 * @param rec {Read} the trimmed record
 * @param opts {Options} the filters to apply
 * @return {bool} true when the record should be kept
 */
export func keeps(rec as fastq.Read, opts as Options) {
    def cfg as Options init normalize($opts);
    if (len($rec.seq) < $cfg.minLength) {
        return false;
    }
    if ($cfg.maxAmbiguous >= 0 and fastq.ambiguousCount($rec) > $cfg.maxAmbiguous) {
        return false;
    }
    if ($cfg.minMeanQuality > 0 and fastq.meanQuality($rec, $cfg.offset) < $cfg.minMeanQuality) {
        return false;
    }
    return true;
}

/**
 * Copy every Nth record of a FASTQ file to a new file.
 *
 * Deterministic and reproducible - the usual way to cut a large run down to a
 * test-sized sample.
 *
 * @param inPath {string} input FASTQ, optionally gzipped
 * @param outPath {string} output FASTQ; gzipped when the name ends in `.gz`
 * @param everyNth {int} keep one read in every N; 1 or less keeps every read
 * @param maxReads {int} stop after writing this many reads; 0 means no limit
 * @return {Stats} counts of what was read and written
 * @throws {Error} when either file cannot be opened, or when `inPath` and
 *   `outPath` are the same file
 */
export func subsampleEvery(inPath as string, outPath as string, everyNth as int, maxReads as int) {
    if (core.sameFile($inPath, $outPath)) {
        core.fail(
            "trim",
            "input and output are the same file (" + $outPath +
                "); the output is opened for writing before the input is read, " +
                "which would truncate it");
    }
    def step as int init $everyNth;
    if ($step < 1) {
        $step = 1;
    }
    def r as fastq.Reader init fastq.open($inPath);
    errdefer quietClose($r);
    def k as source.Sink init source.create($outPath);
    errdefer quietCloseSink($k);
    def seen as int init 0;
    def kept as int init 0;
    def basesIn as int init 0;
    def basesOut as int init 0;
    while (fastq.hasNext($r)) {
        def rec as fastq.Read init fastq.next($r);
        $seen = $seen + 1;
        $basesIn = $basesIn + len($rec.seq);
        if (($seen - 1) % $step != 0) {
            continue;
        }
        if ($maxReads > 0 and $kept >= $maxReads) {
            break;
        }
        fastq.write($k, $rec);
        $kept = $kept + 1;
        $basesOut = $basesOut + len($rec.seq);
    }
    fastq.close($r);
    def code as int init source.closeSink($k);
    if ($code != 0) {
        core.fail(
            "trim",
            "compressor for " + $outPath + " exited with code " + convert.toString($code));
    }
    return countsOnly($seen, $kept, $basesIn, $basesOut);
}

/**
 * Copy a random fraction of a FASTQ file to a new file.
 *
 * Seed the generator with `math.randSeed` first for a reproducible sample.
 *
 * @param inPath {string} input FASTQ, optionally gzipped
 * @param outPath {string} output FASTQ; gzipped when the name ends in `.gz`
 * @param fraction {float} probability of keeping each read, in (0, 1]
 * @return {Stats} counts of what was read and written
 * @throws {Error} when the fraction is out of range, a file cannot be opened, or
 *   `inPath` and `outPath` are the same file
 */
export func subsampleFraction(inPath as string, outPath as string, fraction as float) {
    if ($fraction <= 0.0 or $fraction > 1.0) {
        core.fail("trim", "subsample fraction must be in (0, 1]");
    }
    if (core.sameFile($inPath, $outPath)) {
        core.fail(
            "trim",
            "input and output are the same file (" + $outPath +
                "); the output is opened for writing before the input is read, " +
                "which would truncate it");
    }
    def r as fastq.Reader init fastq.open($inPath);
    errdefer quietClose($r);
    def k as source.Sink init source.create($outPath);
    errdefer quietCloseSink($k);
    def seen as int init 0;
    def kept as int init 0;
    def basesIn as int init 0;
    def basesOut as int init 0;
    while (fastq.hasNext($r)) {
        def rec as fastq.Read init fastq.next($r);
        $seen = $seen + 1;
        $basesIn = $basesIn + len($rec.seq);
        if (math.rand() >= $fraction) {
            continue;
        }
        fastq.write($k, $rec);
        $kept = $kept + 1;
        $basesOut = $basesOut + len($rec.seq);
    }
    fastq.close($r);
    def code as int init source.closeSink($k);
    if ($code != 0) {
        core.fail(
            "trim",
            "compressor for " + $outPath + " exited with code " + convert.toString($code));
    }
    return countsOnly($seen, $kept, $basesIn, $basesOut);
}

/**
 * Fill in the zero value's implicit defaults.
 *
 * @param opts {Options} caller-supplied options
 * @return {Options} the same options with 0 fields resolved to their defaults
 */
func normalize(opts as Options) {
    def cfg as Options init $opts;
    if ($cfg.offset == 0) {
        $cfg.offset = fastq.PHRED33;
    }
    if ($cfg.windowSize < 1) {
        $cfg.windowSize = DEFAULT_WINDOW;
    }
    if ($cfg.minOverlap < 1) {
        $cfg.minOverlap = DEFAULT_OVERLAP;
    }
    return $cfg;
}

/**
 * Record why a trimmed read was discarded, and count it if it was kept.
 *
 * @param st {Stats} the statistics so far
 * @param rec {Read} the trimmed record
 * @param cfg {Options} normalized options
 * @return {Stats} the updated statistics
 */
func judge(st as Stats, rec as fastq.Read, cfg as Options) {
    def out as Stats init $st;
    if (len($rec.seq) < $cfg.minLength) {
        $out.tooShort = $out.tooShort + 1;
        return $out;
    }
    if ($cfg.maxAmbiguous >= 0 and fastq.ambiguousCount($rec) > $cfg.maxAmbiguous) {
        $out.tooAmbiguous = $out.tooAmbiguous + 1;
        return $out;
    }
    if ($cfg.minMeanQuality > 0 and fastq.meanQuality($rec, $cfg.offset) < $cfg.minMeanQuality) {
        $out.lowQuality = $out.lowQuality + 1;
        return $out;
    }
    $out.readsOut = $out.readsOut + 1;
    $out.basesOut = $out.basesOut + len($rec.seq);
    return $out;
}

/**
 * Build a `Stats` for an operation that only copies records.
 *
 * @param seen {int} records read
 * @param kept {int} records written
 * @param basesIn {int} bases read
 * @param basesOut {int} bases written
 * @return {Stats}
 */
func countsOnly(seen as int, kept as int, basesIn as int, basesOut as int) {
    return Stats{
        readsIn: $seen,
        readsOut: $kept,
        basesIn: $basesIn,
        basesOut: $basesOut,
        adapterTrimmed: 0,
        qualityTrimmed: 0,
        tooShort: 0,
        tooAmbiguous: 0,
        lowQuality: 0
    };
}

/**
 * Position of the longest adapter prefix that matches the read's suffix.
 *
 * @param seq {string} the sequence
 * @param probe {string} the adapter
 * @param overlap {int} shortest prefix accepted
 * @return {int} zero-based cut position, or -1 when nothing matches
 */
func partialAt(seq as string, probe as string, overlap as int) {
    def longest as int init len($probe) - 1;
    if ($longest > len($seq)) {
        $longest = len($seq);
    }
    def k as int init $longest;
    while ($k >= $overlap) {
        if (strings.endsWith($seq, strings.substring($probe, 0, $k))) {
            return len($seq) - $k;
        }
        $k = $k - 1;
    }
    return -1;
}

/**
 * Close a reader on a path that is already failing.
 *
 * Cleanup that runs while an error is propagating must not raise one of its own:
 * a deferred throw supersedes the pending error.
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

/**
 * Close a sink on a path that is already failing.
 *
 * @param k {source.Sink} the sink to close
 * @return {int} the close result, or -1 when closing itself failed
 */
func quietCloseSink(k as source.Sink) {
    try {
        return source.closeSink($k);
    } catch (e) {
        return -1;
    }
}
