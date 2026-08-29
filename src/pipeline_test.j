# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * White-box tests for the pipeline orchestration module.
 * Run with: jennifer test pipeline_test.j
 */
use testing;

func nowMillis() {
    def r as os.Result init os.run(["date", "+%s%3N"]);
    return convert.toInt(strings.trim($r.stdout));
}

func testJobRunsAndCaptures() {
    def res as Result init run(job("echo", ["echo", "hello"]));
    testing.assertEqual($res.exitCode, 0);
    testing.assertEqual($res.skipped, false);
    testing.assertEqual(strings.trim($res.stdout), "hello");
}

func testJobArgvIsNotShellInterpreted() {
    def res as Result init run(job("literal", ["echo", "a; echo injected"]));
    testing.assertEqual(strings.trim($res.stdout), "a; echo injected");
}

func testEmptyArgvRejected() {
    testing.assertThrows("emptyArgvJob", "ngs.pipeline");
}

func emptyArgvJob() {
    def none as list of string init [];
    return job("bad", $none);
}

func testJobStdin() {
    def j as Job init feeding(job("cat", ["cat"]), "piped\n");
    def res as Result init run($j);
    testing.assertEqual(strings.trim($res.stdout), "piped");
}

func testShellJobUsesShell() {
    def res as Result init run(shell("pipe", "echo one two | tr ' ' '-'"));
    testing.assertEqual(strings.trim($res.stdout), "one-two");
}

func testFailingJobRaisesOnCheck() {
    testing.assertThrows("failingJob", "ngs.pipeline");
}

func failingJob() {
    def res as Result init run(job("false", ["false"]));
    check($res);
    return $res;
}

func testTolerantJobDoesNotRaise() {
    def res as Result init run(tolerant(job("false", ["false"])));
    check($res);
    testing.assertEqual($res.exitCode, 1);
    testing.assertEqual(allOk([$res]), true);
}

func testAllOkDetectsRealFailure() {
    def res as Result init run(job("false", ["false"]));
    testing.assertEqual(allOk([$res]), false);
}

func testUpToDateJobIsSkippedThenRerun() {
    def input as string init fs.makeTempFile("", "ngsin", ".txt");
    fs.writeString($input, "in\n");
    def output as string init fs.makeTempFile("", "ngsout", ".txt");
    def cmd as string init "cat " + core.shellQuote($input) + " > " + core.shellQuote($output);
    def j as Job init producing(needing(shell("copy", $cmd), [$input]), [$output]);

    def first as Result init run($j);
    testing.assertEqual($first.skipped, true);

    fs.writeString($input, "changed\n");
    def second as Result init run($j);
    testing.assertEqual($second.skipped, false);
    testing.assertEqual($second.exitCode, 0);
    testing.assertEqual(fs.readString($output), "changed\n");

    def third as Result init run($j);
    fs.remove($input);
    fs.remove($output);
    testing.assertEqual($third.skipped, true);
}

func testMissingOutputForcesRun() {
    def input as string init fs.makeTempFile("", "ngsin", ".txt");
    fs.writeString($input, "in\n");
    def output as string init "/tmp/ngs-deck-never-created.txt";
    def j as Job init producing(needing(job("x", ["true"]), [$input]), [$output]);
    def stale as bool init needsRun($j);
    fs.remove($input);
    testing.assertEqual($stale, true);
}

func testJobWithoutOutputsAlwaysRuns() {
    testing.assertEqual(needsRun(job("always", ["true"])), true);
}

func testMissingInputRaises() {
    testing.assertThrows("missingInputJob", "ngs.pipeline");
}

func missingInputJob() {
    def j as Job init producing(
        needing(job("x", ["true"]), ["/nonexistent/ngs/input.txt"]),
        ["/tmp/ngs-deck-never-created.txt"]);
    return needsRun($j);
}

func testRunAllReturnsResultsInOrder() {
    def results as list of Result init runAll([
        job("first", ["echo", "1"]),
        job("second", ["echo", "2"])
    ]);
    testing.assertEqual(len($results), 2);
    testing.assertEqual($results[0].name, "first");
    testing.assertEqual(strings.trim($results[1].stdout), "2");
}

func testRunAllStopsAtFailure() {
    testing.assertThrows("runAllWithFailure", "ngs.pipeline");
}

func runAllWithFailure() {
    return runAll([job("ok", ["true"]), job("bad", ["false"]), job("never", ["true"])]);
}

func testRunParallelKeepsOrderAndRunsConcurrently() {
    def jobs as list of Job init [
        job("one", ["sh", "-c", "sleep 0.4; echo one"]),
        job("two", ["sh", "-c", "sleep 0.4; echo two"]),
        job("three", ["sh", "-c", "sleep 0.4; echo three"]),
        job("four", ["sh", "-c", "sleep 0.4; echo four"])
    ];
    def started as int init nowMillis();
    def results as list of Result init runParallel($jobs, 4);
    def elapsed as int init nowMillis() - $started;
    testing.assertEqual(len($results), 4);
    testing.assertEqual($results[0].name, "one");
    testing.assertEqual($results[3].name, "four");
    testing.assertEqual(strings.trim($results[2].stdout), "three");
    testing.assertEqual(allOk($results), true);
    testing.assertEqual($elapsed < 1200, true);
}

func testRunParallelWidthOneSerialises() {
    def jobs as list of Job init [
        job("a", ["sh", "-c", "sleep 0.3"]),
        job("b", ["sh", "-c", "sleep 0.3"])
    ];
    def started as int init nowMillis();
    def results as list of Result init runParallel($jobs, 1);
    def elapsed as int init nowMillis() - $started;
    testing.assertEqual(len($results), 2);
    testing.assertEqual($elapsed >= 600, true);
}

func testRunParallelHandlesUnevenWaves() {
    def jobs as list of Job init [
        job("a", ["echo", "a"]),
        job("b", ["echo", "b"]),
        job("c", ["echo", "c"]),
        job("d", ["echo", "d"]),
        job("e", ["echo", "e"])
    ];
    def results as list of Result init runParallel($jobs, 2);
    testing.assertEqual(len($results), 5);
    testing.assertEqual($results[4].name, "e");
    testing.assertEqual(strings.trim($results[4].stdout), "e");
}

func testHasToolFindsShell() {
    testing.assertEqual(hasTool("sh"), true);
    testing.assertEqual(hasTool("definitelyNotARealToolXyz"), false);
}

func testHasToolAcceptsPath() {
    testing.assertEqual(hasTool("/nonexistent/binary"), false);
}

func testRequireToolsRaisesOnMissing() {
    testing.assertThrows("requireMissingTool", "ngs.pipeline");
}

func requireMissingTool() {
    return requireTools(["sh", "definitelyNotARealToolXyz"]);
}

func testRequireToolsPassesWhenPresent() {
    requireTools(["sh"]);
    testing.assertEqual(true, true);
}

func testSummaryRendersStatuses() {
    def results as list of Result init runAll([job("ok", ["true"])]);
    testing.assertEqual(strings.contains(summary($results), "ok"), true);
}

func testWorkDirLifecycle() {
    def dir as string init workDir("ngswork");
    testing.assertEqual(fs.isDir($dir), true);
    cleanUp($dir);
    testing.assertEqual(fs.exists($dir), false);
    cleanUp($dir);
}

func testCleanUpRefusesDangerousTargets() {
    testing.assertThrows("cleanRoot", "ngs.pipeline");
    testing.assertThrows("cleanCwd", "ngs.pipeline");
    testing.assertThrows("cleanDot", "ngs.pipeline");
}

func cleanRoot() {
    return cleanUp("/");
}

func cleanCwd() {
    return cleanUp(os.cwd());
}

func cleanDot() {
    return cleanUp(".");
}

func testCleanUpStillRemovesAWorkDir() {
    def dir as string init workDir("ngsguard");
    fs.writeString(path.join($dir, "scratch.txt"), "x");
    cleanUp($dir);
    testing.assertEqual(fs.exists($dir), false);
}
