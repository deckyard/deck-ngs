# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Running external tools, Snakemake-lite. Describe each step as a `Job` with
 * its argv, inputs and outputs; the module skips it when it is up to date,
 * runs it, checks its exit status, and fans independent jobs out with `spawn`.
 *
 * - **Up-to-date checks.** A job whose outputs all exist and are newer than
 *   every input is skipped, so re-running a pipeline resumes it.
 * - **Fan-out with a cap.** `runParallel` runs jobs in waves of at most N,
 *   for when each job is itself a multi-threaded aligner.
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
 * Running jobs needs the `exec` capability (the default `jennifer` binary).
 * Each job's stdout and stderr are captured in memory whole.
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
 * @field name {string} label used in results and `summary`
 * @field argv {list of string} command and arguments, run without a shell
 * @field inputs {list of string} files the job reads
 * @field outputs {list of string} files the job produces
 * @field stdin {string} text fed to standard input, "" for none
 * @field allowFailure {bool} whether a non-zero exit is reported but not fatal
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
 * @field name {string} job label
 * @field exitCode {int} process exit code, 0 when skipped
 * @field stdout {string} captured standard output, "" when skipped
 * @field stderr {string} captured standard error, "" when skipped
 * @field skipped {bool} whether the job was already up to date
 * @field allowFailure {bool} copied from the job
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
 * The command runs without a shell, so arguments need no quoting and cannot be
 * read as shell syntax.
 *
 * @param name {string} label used in results and `summary`
 * @param argv {list of string} command and arguments
 * @return {Job} a job with no inputs, outputs or stdin
 * @throws {Error} kind "ngs.pipeline" when `argv` is empty
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
 * For shell features only: pipes, redirection, globbing. Interpolating
 * untrusted values into `cmdline` is command injection; use `job` with an
 * argv, or quote each value with `core.shellQuote`.
 *
 * @param name {string} label used in results and `summary`
 * @param cmdline {string} POSIX `sh` command line
 * @return {Job} a job that runs `sh -c cmdline`
 */
export func shell(name as string, cmdline as string) {
    return job($name, ["sh", "-c", $cmdline]);
}

/**
 * A copy of the job that declares the files it reads.
 *
 * @param j {Job} the job
 * @param inputs {list of string} files the job reads, replacing any set before
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
 * @param outputs {list of string} files the job writes, replacing any set before
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
 * A job with no declared outputs always runs. Otherwise it runs when an output
 * is missing, or when an input's mtime is at or after the oldest output's: a
 * tie counts as stale, since coarse filesystem timestamps make ties common.
 *
 * @param j {Job} the job
 * @return {bool} true when the job is out of date
 * @throws {Error} kind "ngs.pipeline" when a declared input does not exist and
 *   the job has outputs
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
 * Run one job, skipping it when `needsRun` says it is up to date. A non-zero
 * exit is returned, not raised; pass the result to `check`.
 *
 * @param j {Job} job to run
 * @return {Result} what happened
 * @throws {Error} kind "ngs.pipeline" without the `exec` capability or when an
 *   input is missing
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
 * @param jobs {list of Job} jobs in dependency order
 * @return {list of Result} one result per job, in order
 * @throws {Error} kind "ngs.pipeline" on the first job that fails without
 *   `allowFailure`; the results so far are lost
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
 * Jobs run in waves: a wave finishes before the next begins. For independent
 * per-sample work; it does not resolve dependencies between jobs. Failures are
 * returned, not raised; pass the results to `checkAll`.
 *
 * @param jobs {list of Job} independent jobs
 * @param width {int} most jobs running at once; 0 or less for `os.NCPU`
 * @return {list of Result} results in the order of `jobs`
 * @throws {Error} kind "ngs.pipeline" without the `exec` capability or when a
 *   job's input is missing
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
 * @param res {Result} result to inspect
 * @throws {Error} kind "ngs.pipeline" when the job exited non-zero without
 *   `allowFailure`; the message carries the job's stderr
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
 * @throws {Error} kind "ngs.pipeline" on the first job that failed without
 *   `allowFailure`
 */
export func checkAll(results as list of Result) {
    for (def res in $results) {
        check($res);
    }
}

/**
 * Whether every result succeeded or was tolerated: `checkAll` without raising.
 *
 * @param results {list of Result} results to inspect
 * @return {bool} true when no real failure is present
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
 * Whether a program is on `PATH`; see `core.hasTool`. Works without `exec`,
 * and checks only that a regular file exists, not that it is executable.
 *
 * @param name {string} program name, or a path containing a separator
 * @return {bool} true when a regular file of that name is found
 */
export func hasTool(name as string) {
    return core.hasTool($name);
}

/**
 * Raise unless every named tool is available.
 *
 * Call it at the top of a pipeline so a missing dependency fails before the
 * first job, not halfway through.
 *
 * @param names {list of string} program names the pipeline needs
 * @throws {Error} kind "ngs.pipeline" listing every missing tool
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
 * A fresh scratch directory under the system temp directory. The caller
 * removes it, usually with `defer pipeline.cleanUp(dir)`.
 *
 * @param tag {string} directory name prefix
 * @return {string} path of the new directory
 */
export func workDir(tag as string) {
    return fs.makeTempDir("", $tag);
}

/**
 * Remove a directory and everything under it.
 *
 * Pair with `workDir` under a `defer`. Refuses "", `.`, `..`, the filesystem
 * root and the working directory, so a mistyped variable cannot become a
 * recursive delete of something that matters. Anything else is removed
 * without further checks; a path that is not a directory is left alone.
 *
 * @param dir {string} directory to remove
 * @throws {Error} kind "ngs.pipeline" for a refused target
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
 * Render results as one `<status> <name>` line per job; status is `skip`,
 * `ok` (padded to four characters), `warn` for a tolerated failure, or
 * `fail`.
 *
 * @param results {list of Result} results to report
 * @return {string} the lines, each newline-terminated
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
 * @param j {Job} job to execute
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
 * @return {string} "skip", "ok  " (padded to four characters), "warn" or "fail"
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
