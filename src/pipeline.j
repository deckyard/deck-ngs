# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Pipeline orchestration for NGS tools - a Snakemake-lite in a real language.
 *
 * NGS analysis is fundamentally `bwa | samtools | gatk` glued together. This
 * module is the glue: describe each step as a `Job` with its inputs, outputs and
 * argv, and the module decides whether it needs to run, runs it, checks its exit
 * status, and fans independent jobs out across cores with `spawn`.
 *
 * Two ideas do most of the work:
 *
 * - **Up-to-date checks.** A job whose outputs all exist and are newer than
 *   every input is skipped, so re-running a pipeline resumes rather than redoes.
 * - **Fan-out with a concurrency cap.** `runParallel` runs jobs in waves of at
 *   most N, which is what you want when each job is itself a multi-threaded
 *   aligner.
 *
 * ```jennifer
 * def jobs as list of pipeline.Job init [];
 * for (def s in $samples) {
 *     def j as pipeline.Job init pipeline.job($s, ["bwa", "mem", "ref.fa", $s + ".fq"]);
 *     $jobs[] = pipeline.producing(pipeline.needing($j, [$s + ".fq"]), [$s + ".bam"]);
 * }
 * def results as list of pipeline.Result init pipeline.runParallel($jobs, 4);
 * pipeline.checkAll($results);
 * ```
 *
 * Everything here needs the `exec` capability, so the default `jennifer` binary.
 *
 * @module pipeline
 * @since 0.1.0
 */
use os;
use fs;
use io;
use task;
use meta;
use path;
use strings;
use convert;
import "./core.j" as core;

/**
 * One step of a pipeline.
 *
 * @field name {string} label used in logs and results
 * @field argv {list of string} the command and its arguments, run without a shell
 * @field inputs {list of string} files this job reads
 * @field outputs {list of string} files this job produces
 * @field stdin {string} text fed to the process's standard input
 * @field allowFailure {bool} when true, a non-zero exit is reported but not fatal
 */
export def struct Job {
    name as string,
    argv as list of string,
    inputs as list of string,
    outputs as list of string,
    stdin as string,
    allowFailure as bool
};

/**
 * What running a job produced.
 *
 * @field name {string} the job's label
 * @field exitCode {int} the process exit code, 0 when skipped
 * @field stdout {string} captured standard output
 * @field stderr {string} captured standard error
 * @field skipped {bool} true when the job was already up to date
 * @field allowFailure {bool} carried over from the job
 */
export def struct Result {
    name as string,
    exitCode as int,
    stdout as string,
    stderr as string,
    skipped as bool,
    allowFailure as bool
};

/**
 * Build a job from a name and an argv.
 *
 * The command is executed directly, not through a shell, so arguments never need
 * quoting and cannot be reinterpreted as shell syntax.
 *
 * @param name {string} label used in logs and results
 * @param argv {list of string} the command and its arguments
 * @return {Job} a job with no declared inputs or outputs
 * @throws {Error} when argv is empty
 */
export func job(name as string, argv as list of string) {
    if (len($argv) == 0) {
        core.fail("pipeline", "job " + $name + " has an empty argv");
    }
    def none as list of string init [];
    return Job{
        name: $name,
        argv: $argv,
        inputs: $none,
        outputs: $none,
        stdin: "",
        allowFailure: false
    };
}

/**
 * Build a job that runs a shell command line.
 *
 * Use this only when you genuinely need shell features - pipes, redirection,
 * globbing. Interpolating untrusted values into `cmdline` is a command-injection
 * risk; `job` with an argv is the safe default, and `core.shellQuote` is
 * available when you must build a command line by hand.
 *
 * @param name {string} label used in logs and results
 * @param cmdline {string} a POSIX `sh` command line
 * @return {Job} a job that runs `sh -c cmdline`
 */
export func shell(name as string, cmdline as string) {
    return job($name, ["sh", "-c", $cmdline]);
}

/**
 * A copy of the job that declares the files it reads.
 *
 * @param j {Job} the job
 * @param inputs {list of string} files the job reads
 * @return {Job} a copy with `inputs` set
 */
export func needing(j as Job, inputs as list of string) {
    def out as Job init $j;
    $out.inputs = $inputs;
    return $out;
}

/**
 * A copy of the job that declares the files it produces.
 *
 * @param j {Job} the job
 * @param outputs {list of string} files the job writes
 * @return {Job} a copy with `outputs` set
 */
export func producing(j as Job, outputs as list of string) {
    def out as Job init $j;
    $out.outputs = $outputs;
    return $out;
}

/**
 * A copy of the job with text to feed to its standard input.
 *
 * @param j {Job} the job
 * @param text {string} data written to the process's stdin
 * @return {Job} a copy with `stdin` set
 */
export func feeding(j as Job, text as string) {
    def out as Job init $j;
    $out.stdin = $text;
    return $out;
}

/**
 * A copy of the job whose failure is reported but not fatal.
 *
 * @param j {Job} the job
 * @return {Job} a copy with `allowFailure` set
 */
export func tolerant(j as Job) {
    def out as Job init $j;
    $out.allowFailure = true;
    return $out;
}

/**
 * Whether a job still needs to run.
 *
 * A job with no declared outputs always runs. Otherwise it runs when any output
 * is missing, or when any input is newer than the oldest output.
 *
 * @param j {Job} the job
 * @return {bool} true when the job is out of date
 * @throws {Error} when a declared input does not exist
 */
export func needsRun(j as Job) {
    if (len($j.outputs) == 0) {
        return true;
    }
    for (def want in $j.inputs) {
        if (not fs.exists($want)) {
            core.fail("pipeline", "job " + $j.name + " needs missing input: " + $want);
        }
    }
    def oldest as int init -1;
    for (def made in $j.outputs) {
        if (not fs.exists($made)) {
            return true;
        }
        def stamp as int init fs.stat($made).mtimeNanos;
        if ($oldest < 0 or $stamp < $oldest) {
            $oldest = $stamp;
        }
    }
    for (def want in $j.inputs) {
        # `>=`, not `>`. An input stamped the same as the oldest output cannot
        # be shown to predate it, so the output cannot be shown to have been
        # built from it - and on a filesystem whose timestamps are coarser than
        # the time a job takes, that is the common case rather than a rare tie.
        # Re-running a job that was in fact current costs time; skipping one
        # that was stale ships a wrong result quietly, so the tie breaks
        # towards running.
        if (fs.stat($want).mtimeNanos >= $oldest) {
            return true;
        }
    }
    return false;
}

/**
 * Run one job, skipping it when its outputs are already up to date.
 *
 * @param j {Job} the job to run
 * @return {Result} what happened
 * @throws {Error} without the `exec` capability, or when an input is missing
 */
export func run(j as Job) {
    requireExec();
    if (not needsRun($j)) {
        return Result{
            name: $j.name,
            exitCode: 0,
            stdout: "",
            stderr: "",
            skipped: true,
            allowFailure: $j.allowFailure
        };
    }
    def r as os.Result init execute($j);
    return Result{
        name: $j.name,
        exitCode: $r.exitCode,
        stdout: $r.stdout,
        stderr: $r.stderr,
        skipped: false,
        allowFailure: $j.allowFailure
    };
}

/**
 * Run jobs one after another, stopping at the first real failure.
 *
 * @param jobs {list of Job} the jobs, in dependency order
 * @return {list of Result} results for the jobs that ran, in order
 * @throws {Error} on the first job that fails without `allowFailure`
 */
export func runAll(jobs as list of Job) {
    def out as list of Result init [];
    for (def j in $jobs) {
        def res as Result init run($j);
        $out[] = $res;
        check($res);
    }
    return $out;
}

/**
 * Run jobs concurrently, at most `width` at a time.
 *
 * Jobs are launched in waves, so a wave finishes before the next begins. Use
 * this for independent per-sample work; it does not resolve dependencies between
 * the jobs you hand it.
 *
 * @param jobs {list of Job} independent jobs
 * @param width {int} maximum jobs running at once; 0 or less means `os.NCPU`
 * @return {list of Result} results in the same order as `jobs`
 * @throws {Error} without the `exec` capability
 */
export func runParallel(jobs as list of Job, width as int) {
    requireExec();
    def lanes as int init $width;
    if ($lanes < 1) {
        $lanes = os.NCPU;
    }
    def out as list of Result init [];
    def at as int init 0;
    while ($at < len($jobs)) {
        def stop as int init $at + $lanes;
        if ($stop > len($jobs)) {
            $stop = len($jobs);
        }
        def wave as list of task of Result init [];
        for (def i in $at..$stop) {
            def one as Job init $jobs[$i];
            $wave[] = spawn {
                return run($one);
            };
        }
        for (def res in task.waitAll($wave)) {
            $out[] = $res;
        }
        $at = $stop;
    }
    return $out;
}

/**
 * Raise when a result represents a real failure.
 *
 * @param res {Result} the result to inspect
 * @throws {Error} when the job exited non-zero without `allowFailure`
 */
export func check(res as Result) {
    if ($res.exitCode == 0 or $res.allowFailure) {
        return;
    }
    def code as string init convert.toString($res.exitCode);
    core.fail("pipeline", "job " + $res.name + " failed with code " + $code + ": " + $res.stderr);
}

/**
 * Raise on the first real failure in a batch of results.
 *
 * @param results {list of Result} results to inspect
 * @throws {Error} on the first job that failed without `allowFailure`
 */
export func checkAll(results as list of Result) {
    for (def res in $results) {
        check($res);
    }
}

/**
 * Whether every result succeeded or was tolerated.
 *
 * @param results {list of Result} results to inspect
 * @return {bool}
 */
export func allOk(results as list of Result) {
    for (def res in $results) {
        if ($res.exitCode != 0 and not $res.allowFailure) {
            return false;
        }
    }
    return true;
}

/**
 * Whether an executable is on `PATH`.
 *
 * Walks `PATH` with `fs` rather than shelling out, so it works without the
 * `exec` capability and cannot be confused by shell quoting.
 *
 * @param name {string} program name, or an absolute path
 * @return {bool} true when the program exists and is a file
 */
export func hasTool(name as string) {
    return core.hasTool($name);
}

/**
 * Raise unless every named tool is available.
 *
 * Call this once at the top of a pipeline so a missing dependency fails
 * immediately with a clear message, rather than halfway through a long run.
 *
 * @param names {list of string} program names the pipeline needs
 * @throws {Error} listing every tool that is missing
 */
export func requireTools(names as list of string) {
    def missing as list of string init [];
    for (def name in $names) {
        if (not hasTool($name)) {
            $missing[] = $name;
        }
    }
    if (len($missing) > 0) {
        core.fail("pipeline", "missing required tools on PATH: " + strings.join($missing, ", "));
    }
}

/**
 * A scratch directory that the caller is responsible for removing.
 *
 * @param tag {string} prefix used in the directory name
 * @return {string} path to a fresh directory
 */
export func workDir(tag as string) {
    return fs.makeTempDir("", $tag);
}

/**
 * Remove a directory and everything under it.
 *
 * Pair with `workDir` under a `defer` so intermediates are cleaned up on every
 * exit path.
 *
 * Refuses a filesystem root, the working directory and the relative `.` / `..`,
 * so a mistyped variable cannot turn a cleanup into a recursive delete of
 * something that matters.
 *
 * @param dir {string} directory to remove
 * @throws {Error} when asked to remove a root or the working directory
 */
export func cleanUp(dir as string) {
    def target as string init path.clean($dir);
    if (len($target) == 0 or $target == "." or $target == ".." or $target == os.DIRSEP) {
        core.fail("pipeline", "refusing to recursively remove " + $target);
    }
    if ($target == path.clean(os.cwd())) {
        core.fail("pipeline", "refusing to recursively remove the working directory");
    }
    if (fs.isDir($target)) {
        fs.removeAll($target);
    }
}

/**
 * Render results as one status line per job.
 *
 * @param results {list of Result} results to report
 * @return {string} a multi-line summary, newline-terminated
 */
export func summary(results as list of Result) {
    def out as string init "";
    for (def res in $results) {
        $out = $out + io.sprintf("%s %s\n", statusOf($res), $res.name);
    }
    return $out;
}

/**
 * Run a job's process, with or without stdin.
 *
 * @param j {Job} the job to execute
 * @return {os.Result} the captured process result
 */
func execute(j as Job) {
    if (len($j.stdin) == 0) {
        return os.run($j.argv);
    }
    return os.run($j.argv, $j.stdin);
}

/**
 * The one-word status of a result.
 *
 * @param res {Result} the result
 * @return {string} "skip", "ok", "warn" or "fail"
 */
func statusOf(res as Result) {
    if ($res.skipped) {
        return "skip";
    }
    if ($res.exitCode == 0) {
        return "ok  ";
    }
    if ($res.allowFailure) {
        return "warn";
    }
    return "fail";
}

/**
 * Raise when the interpreter cannot start processes.
 *
 * @throws {Error} without the `exec` capability
 */
func requireExec() {
    if (not meta.hasCapability("exec")) {
        core.fail(
            "pipeline",
            "running external tools needs the `exec` capability (default jennifer binary)");
    }
}
