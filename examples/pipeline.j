#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0
#
# A Snakemake-lite per-sample pipeline: QC each sample, trim it, then run an
# external aligner - fanned out across cores, resumable, and failing loudly.
#
#   jennifer run examples/pipeline.j sampleA.fastq.gz sampleB.fastq.gz
#
# The alignment step here shells out to `wc -c` as a stand-in so the example runs
# anywhere; swap it for the real `bwa mem` argv and the structure is unchanged.
# That substitution is the point of the deck: alignment stays with the native
# tool, and this file is the glue around it.
use io;
use os;
use fs;
use path;
use strings;
import "../src/qc.j" as qc;
import "../src/trim.j" as trim;
import "../src/pipeline.j" as pipeline;
import "../src/core.j" as core;

if (len(os.ARGS) < 2) {
    io.eprintf("usage: pipeline.j <sample.fastq[.gz]> [more samples...]\n");
    exit 2;
}

# Fail immediately if the tools this pipeline drives are not installed, rather
# than halfway through a long run.
pipeline.requireTools(["gzip", "wc"]);

def work as string init pipeline.workDir("ngspipe");
defer pipeline.cleanUp($work);
io.printf("work directory %s\n", $work);

def samples as list of string init [];
for (def i in 1..len(os.ARGS)) {
    $samples[] = os.ARGS[$i];
}

# --- stage 1: QC and trim each sample, in this interpreter --------------------

def trimmed as list of string init [];
for (def sample in $samples) {
    def label as string init path.stem(path.base($sample));
    def opts as qc.Options init qc.options();
    $opts.maxReads = 100000;
    def report as qc.Report init qc.analyze($sample, $opts);
    io.printf("\n== %s ==\n", $label);
    io.printf("%s", qc.summary($report));

    def out as string init path.join($work, $label + ".trimmed.fastq.gz");
    def st as trim.Stats init trim.run($sample, $out, trim.options());
    io.printf("trimmed %d of %d reads -> %s\n", $st.readsOut, $st.readsIn, $out);
    $trimmed[] = $out;
}

# --- stage 2: fan the external step out across cores --------------------------

def jobs as list of pipeline.Job init [];
for (def input in $trimmed) {
    def label as string init path.base($input);
    def target as string init $input + ".count";
    def cmd as string init "gzip -dc " + quoted($input) + " | wc -l > " + quoted($target);
    def j as pipeline.Job init pipeline.shell($label, $cmd);
    $jobs[] = pipeline.producing(pipeline.needing($j, [$input]), [$target]);
}

io.printf("\nrunning %d jobs across %d cores\n", len($jobs), os.NCPU);
def results as list of pipeline.Result init pipeline.runParallel($jobs, os.NCPU);
io.printf("%s", pipeline.summary($results));
pipeline.checkAll($results);

# Re-running is a no-op: every output is newer than its input.
def again as list of pipeline.Result init pipeline.runParallel($jobs, os.NCPU);
io.printf("second pass (should all skip)\n%s", pipeline.summary($again));

for (def input in $trimmed) {
    def counted as string init fs.readString($input + ".count");
    io.printf("%s -> %s lines\n", path.base($input), trimText($counted));
}

func quoted(p as string) {
    return core.shellQuote($p);
}

func trimText(s as string) {
    return strings.trim($s);
}
