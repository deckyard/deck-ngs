# Pipeline orchestration

This is the deck's sweet spot. NGS analysis is fundamentally `bwa | samtools |
gatk` glued together, and Jennifer already has `os.run`, `os.spawn` and `spawn`
concurrency. The `pipeline` module is the glue: describe each step as a `Job`
with its inputs, outputs and argv, and it decides whether the step needs to run,
runs it, checks its exit status, and fans independent steps out across cores.

A Snakemake-lite, in a real language.

```jennifer
import "@jennifer/ngs/pipeline.j" as pipeline;

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

Everything here needs the `exec` capability — the default `jennifer` binary.

## Building a job

`pipeline.job(name, argv)` runs a command **without a shell**. Arguments never
need quoting and cannot be reinterpreted as shell syntax, which is what you want
for anything built from a filename or a sample identifier.

The builders are copy-returning, so they chain in any order:

```jennifer
def j as pipeline.Job init pipeline.job("align", ["bwa", "mem", "ref.fa", "r.fq"]);
$j = pipeline.needing($j, ["ref.fa", "r.fq"]);   # files it reads
$j = pipeline.producing($j, ["out.bam"]);        # files it writes
$j = pipeline.feeding($j, $text);                # data for its stdin
$j = pipeline.tolerant($j);                      # a non-zero exit is not fatal
```

`pipeline.shell(name, cmdline)` is there when you genuinely need shell features
— a pipe, a redirect, a glob:

```jennifer
def cmd as string init "bwa mem ref.fa " + core.shellQuote($fq)
    + " | samtools sort -o " + core.shellQuote($bam);
def j as pipeline.Job init pipeline.shell($sample, $cmd);
```

Interpolating an untrusted value into `cmdline` is a command-injection risk.
`core.shellQuote` is the tool for building one safely; `pipeline.job` with an
argv is the safe default.

## Skipping work that is already done

A job that declares outputs is **skipped when every output exists and is newer
than every input**. Re-running a pipeline therefore resumes rather than redoes,
which is the property that makes a long pipeline survivable.

```jennifer
pipeline.needsRun($j);      # the decision, on its own
def res as pipeline.Result init pipeline.run($j);
$res.skipped;               # true when nothing ran
```

A job with **no declared outputs always runs** — with nothing to compare, there
is no basis for skipping.

A declared input that does not exist raises `ngs.pipeline` rather than silently
running a step that cannot work.

## Running

```jennifer
pipeline.run($j);                    # one job
pipeline.runAll($jobs);              # in order, stopping at the first failure
pipeline.runParallel($jobs, 8);      # concurrently, at most 8 at a time
```

`runParallel` launches jobs in **waves**: a wave finishes before the next
begins. Results come back in the same order as the input, regardless of the
order they finished. Pass `0` for the width to use `os.NCPU`.

It runs independent jobs. It does not resolve dependencies between the jobs you
hand it — order those yourself, or run them in stages:

```jennifer
def aligned as list of pipeline.Result init pipeline.runParallel($alignJobs, 4);
pipeline.checkAll($aligned);
def called as list of pipeline.Result init pipeline.runParallel($callJobs, 4);
```

Width is worth thinking about: each job here is usually itself a multi-threaded
aligner. Four jobs of four threads each will beat sixteen jobs of one.

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
pipeline.summary($results);    # one "ok / skip / warn / fail" line per job
```

`runAll` checks as it goes, so it stops at the first failure. `runParallel` does
not — a wave is already running — so check the results yourself.

Note that stdout and stderr are **captured into memory**. That is fine for tool
logs, which is what they are. For a tool that writes a large result to stdout,
redirect it to a file with `pipeline.shell` rather than capturing it.

## Requirements and scratch space

Check for tools once, at the top, so a missing dependency fails immediately
rather than after an hour of alignment:

```jennifer
pipeline.requireTools(["bwa", "samtools", "gatk"]);
pipeline.hasTool("bwa");     # the same test, as a bool
```

`hasTool` walks `PATH` with `fs` rather than shelling out, so it works without
the `exec` capability and cannot be confused by quoting.

For intermediates, pair `workDir` with `defer` so they are cleaned up on every
exit path — including a `throw`:

```jennifer
def work as string init pipeline.workDir("myrun");
defer pipeline.cleanUp($work);
```

## A complete example

`examples/pipeline.j` in the repository is a runnable version of the shape
above: it QCs and trims each sample in-process, then fans an external step out
across cores, then demonstrates that a second pass skips everything. See
[Recipes](recipes.md#a-per-sample-pipeline) for the walk-through.
