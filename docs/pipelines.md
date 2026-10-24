# Pipeline orchestration

The `pipeline` module runs external tools in the manner of a small Snakemake.
Each step is a `Job` with an argv, its inputs and its outputs; the module
decides whether the step needs to run, runs it, checks its exit status and fans
independent steps out across cores with `spawn`.

```jennifer
import "@deckyard/ngs/pipeline.j" as pipeline;

pipeline.requireTools(["bwa", "samtools"]);

def jobs as list of pipeline.Job init [];
for (def s in $samples) {
    def j as pipeline.Job init pipeline.job($s, ["bwa", "mem", "ref.fa", $s + ".fq"]);
    $jobs[] = pipeline.producing(pipeline.needing($j, [$s + ".fq"]), [$s + ".bam"]);
}

def results as list of pipeline.Result init pipeline.runParallel($jobs, 8);
pipeline.checkAll($results);
io.printf("%s", pipeline.summary($results));
```

`pipeline.run`, `runAll` and `runParallel` need the `exec` capability (the
default `jennifer` binary) and raise `ngs.pipeline` without it. Every error
here is `ngs.pipeline`.

## Building a job

`pipeline.job(name, argv)` runs the command without a shell, so arguments need
no quoting and cannot be read as shell syntax. An empty argv raises
`ngs.pipeline`.

The builders return a modified copy, so they chain in any order:

```jennifer
def j as pipeline.Job init pipeline.job("align", ["bwa", "mem", "ref.fa", "r.fq"]);
$j = pipeline.needing($j, ["ref.fa", "r.fq"]);   # files it reads
$j = pipeline.producing($j, ["out.bam"]);        # files it writes
$j = pipeline.feeding($j, $text);                # text for its stdin
$j = pipeline.tolerant($j);                      # a non-zero exit is not fatal
```

`pipeline.shell(name, cmdline)` runs `sh -c cmdline`, for pipes, redirects and
globs:

```jennifer
def cmd as string init "bwa mem ref.fa " + core.shellQuote($fq)
    + " | samtools sort -o " + core.shellQuote($bam);
def j as pipeline.Job init pipeline.shell($sample, $cmd);
```

Interpolating an untrusted value into `cmdline` is a command-injection risk.
Quote every value with `core.shellQuote`, or use `pipeline.job` with an argv.

## Skipping finished work

A job that declares outputs is skipped when every output exists and is strictly
newer than every input, so re-running a pipeline resumes it.

An input with the same timestamp as the oldest output counts as stale: equal
timestamps cannot show which was written first, and on a filesystem whose
timestamp resolution is coarser than a job's run time a whole pipeline can land
in one tick. Re-running a current job costs time; skipping a stale one ships a
wrong result without an error.

```jennifer
pipeline.needsRun($j);      # the decision on its own
def res as pipeline.Result init pipeline.run($j);
$res.skipped;               # true when nothing ran; exitCode is then 0
```

A job with no declared outputs always runs. A declared input that does not
exist raises `ngs.pipeline`.

## Running

```jennifer
pipeline.run($j);                    # one job
pipeline.runAll($jobs);              # in order, stopping at the first failure
pipeline.runParallel($jobs, 8);      # concurrently, at most 8 at a time
```

`runParallel` launches jobs in waves of `width`; a wave finishes before the
next starts. Results come back in input order. A width of 0 or less means
`os.NCPU`.

`runParallel` does not resolve dependencies between the jobs it is given. Run
dependent steps as separate stages:

```jennifer
def aligned as list of pipeline.Result init pipeline.runParallel($alignJobs, 4);
pipeline.checkAll($aligned);
def called as list of pipeline.Result init pipeline.runParallel($callJobs, 4);
```

Each job is often a multi-threaded aligner itself: four jobs of four threads
each outperform sixteen jobs of one.

## Results and failure

```jennifer
def struct Result {
    name as string,
    exitCode as int,
    stdout as string,
    stderr as string,
    skipped as bool,
    allowFailure as bool
};
```

```jennifer
pipeline.check($res);          # raises unless the job succeeded or is tolerant
pipeline.checkAll($results);   # raises on the first real failure
pipeline.allOk($results);      # the same question, as a bool
pipeline.summary($results);    # one "skip", "ok", "warn" or "fail" line per job
```

The `ngs.pipeline` error from `check` carries the job's name, exit code and
stderr. `runAll` checks as it goes and stops at the first failure.
`runParallel` does not check, since the wave is already running; check its
results yourself.

stdout and stderr are captured into memory. For a tool that writes a large
result to stdout, redirect it to a file with `pipeline.shell`.

## Requirements and scratch space

Check for tools at the top, so a missing one fails before an hour of
alignment:

```jennifer
pipeline.requireTools(["bwa", "samtools", "gatk"]);   # raises, listing every missing tool
pipeline.hasTool("bwa");                              # the same test, as a bool
```

`hasTool` walks `PATH` with `fs` rather than running a shell, so it works
without `exec`. A name containing a path separator is checked as a path.

Pair `workDir` with `defer` so intermediates are removed on every exit path,
including a `throw`:

```jennifer
def work as string init pipeline.workDir("myrun");
defer pipeline.cleanUp($work);
```

`cleanUp` raises `ngs.pipeline` rather than recursively delete a filesystem
root, the working directory, `.` or `..`.

`examples/pipeline.j` is a runnable version: it QCs and trims each sample
in-process, fans an external step out across cores, and shows a second pass
skipping everything. [Recipes](recipes.md#a-per-sample-pipeline) walks through
it.
