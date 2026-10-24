# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming line sources and sinks. Every reader in the deck pulls one line
 * at a time from a `Source` and every writer pushes one line at a time into a
 * `Sink`, so a 50 GB FASTQ never lands in memory.
 *
 * Plain files stream through an `fs.File` handle. Gzipped files stream through
 * an external `gzip` whose output goes into a FIFO opened as an ordinary read
 * handle: `compress` has a streaming compressor but no streaming
 * decompressor. Memory stays at one line whatever the file size; the cost is
 * the `exec` capability (the default `jennifer` binary), a Unix host, and
 * `gzip` and `mkfifo` on `PATH`.
 *
 * `openGzipBuffered` is the fallback for hosts without `exec`. It holds the
 * decompressed file in memory, so use it for small inputs only.
 *
 * @module source
 * @since 0.1.0
 */
use fs;
use os;
use meta;
use path;
use strings;
use compress;
use convert;
import "./core.j" as core;

/** First byte of the gzip magic number (RFC 1952). */
export def const GZIP_MAGIC1 as int init 0x1f;

/** Second byte of the gzip magic number (RFC 1952). */
export def const GZIP_MAGIC2 as int init 0x8b;

/**
 * An open line source. Treat it as an opaque handle: the `fs.File` and
 * `os.Process` fields are registry handles, so copies share one stream.
 *
 * @field file {fs.File} readable handle lines are pulled from
 * @field proc {os.Process} decompressor child, when `piped` is true
 * @field helper {string} private directory or temp file removed on close, "" when none
 * @field piped {bool} whether a child process feeds this source
 * @field origin {string} path opened, or "<stdin>", for messages
 */
export def struct Source {
    file as fs.File,
    proc as os.Process,
    helper as string,
    piped as bool,
    origin as string
};

/**
 * An open line sink, the write side of `Source`.
 *
 * @field file {fs.File} writable handle lines are pushed into
 * @field proc {os.Process} compressor child, when `piped` is true
 * @field helper {string} private directory removed on close, "" when none
 * @field piped {bool} whether a child process drains this sink
 * @field origin {string} path written, for messages
 */
export def struct Sink {
    file as fs.File,
    proc as os.Process,
    helper as string,
    piped as bool,
    origin as string
};

/**
 * Whether a file starts with the gzip magic number.
 *
 * Sniffs the content, not the extension, so a gzipped `.fq` and a plain `.gz`
 * are both read correctly.
 *
 * @param path {string} file to inspect
 * @return {bool} true when the first two bytes are 1f 8b
 * @throws {Error} kind "ngs.source" when the file does not exist
 */
export func isGzip(path as string) {
    if (not fs.exists($path)) {
        core.fail("source", "no such file: " + $path);
    }
    def f as fs.File init fs.open($path, "read");
    def head as bytes init fs.readBytes($f, 2);
    fs.close($f);
    if (len($head) < 2) {
        return false;
    }
    return $head[0] == GZIP_MAGIC1 and $head[1] == GZIP_MAGIC2;
}

/**
 * Open a file for streaming line reads, decompressing gzip transparently.
 *
 * @param path {string} file to read
 * @return {Source} an open source; close it with `close` or `closeStrict`
 * @throws {Error} kind "ngs.source" when the file is missing, or is gzipped and
 *   the host lacks `exec` or `mkfifo` fails
 */
export func open(path as string) {
    if (isGzip($path)) {
        return openGzip($path);
    }
    return openPlain($path);
}

/**
 * Open an uncompressed file for streaming line reads.
 *
 * @param path {string} file to read
 * @return {Source} an open source
 * @throws {Error} kind "ngs.source" when the file is missing
 */
export func openPlain(path as string) {
    if (not fs.exists($path)) {
        core.fail("source", "no such file: " + $path);
    }
    def f as fs.File init fs.open($path, "read");
    def p as os.Process;
    return Source{file: $f, proc: $p, helper: "", piped: false, origin: $path};
}

/**
 * Open standard input as a line source, for Unix pipes.
 *
 * Reopens `/dev/stdin`, so it is Unix-only, and only works when stdin can be
 * reopened: a pipe, a file or a terminal, not a socket.
 *
 * @return {Source} an open source over the process's stdin, `origin` "<stdin>"
 * @throws {Error} kind "ngs.source" when `/dev/stdin` does not exist; kind
 *   "runtime" when it exists but cannot be opened (a socket gives ENXIO)
 */
export func openStdin() {
    if (not fs.exists("/dev/stdin")) {
        core.fail("source", "stdin streaming needs /dev/stdin (Unix hosts only)");
    }
    def f as fs.File init fs.open("/dev/stdin", "read");
    def p as os.Process;
    return Source{file: $f, proc: $p, helper: "", piped: false, origin: "<stdin>"};
}

/**
 * Open a gzipped file for streaming line reads at constant memory.
 *
 * Spawns `sh -c 'exec gzip -dc -- <path> > <fifo>'` with the FIFO in a private
 * directory, then reads the FIFO as a normal line handle. The path is
 * shell-quoted. `close` returns gzip's exit code.
 *
 * @param path {string} gzipped file to read
 * @return {Source} an open source
 * @throws {Error} kind "ngs.source" when the file is missing, without the `exec`
 *   capability, or when `mkfifo` fails
 */
export func openGzip(path as string) {
    if (not fs.exists($path)) {
        core.fail("source", "no such file: " + $path);
    }
    if (not meta.hasCapability("exec")) {
        core.fail(
            "source",
            "streaming gzip needs the `exec` capability (default jennifer binary); " +
                "use source.openGzipBuffered for a small-file fallback");
    }
    def work as string init makeFifoDir("ngsread");
    errdefer discard($work);
    def fifo as string init fifoIn($work);
    def cmd as string init "exec gzip -dc -- " + core.shellQuote($path) + " > " +
        core.shellQuote($fifo);
    def child as os.Process init os.spawn(["sh", "-c", $cmd]);
    errdefer reap($child);
    def f as fs.File init fs.open($fifo, "read");
    return Source{file: $f, proc: $child, helper: $work, piped: true, origin: $path};
}

/**
 * Open a gzipped file by decompressing it whole, for hosts without `exec`.
 *
 * Holds the whole decompressed file in memory, then writes it to a temp file
 * that `close` removes. For small inputs only: a small crafted `.gz` can
 * expand to gigabytes. `openGzip` is the streaming path.
 *
 * @param path {string} gzipped file to read
 * @return {Source} an open source over a temporary decompressed copy
 * @throws {Error} kind "ngs.source" when the file is missing; the `compress`
 *   error when it is not valid gzip
 */
export func openGzipBuffered(path as string) {
    if (not fs.exists($path)) {
        core.fail("source", "no such file: " + $path);
    }
    def packed as bytes init fs.readBytes($path);
    def raw as bytes init compress.unpack($packed, "gzip");
    def temp as string init fs.makeTempFile("", "ngsbuf", ".txt");
    errdefer discard($temp);
    fs.writeBytes($temp, $raw);
    def f as fs.File init fs.open($temp, "read");
    def p as os.Process;
    return Source{file: $f, proc: $p, helper: $temp, piped: false, origin: $path};
}

/**
 * Whether another line can be read.
 *
 * @param s {Source} an open source
 * @return {bool} false once the stream is exhausted
 */
export func hasNext(s as Source) {
    return not fs.eof($s.file);
}

/**
 * Read the next line, with its line terminator stripped.
 *
 * @param s {Source} an open source
 * @return {string} the next line
 * @throws {Error} kind "runtime" when the source is exhausted; guard with `hasNext`
 */
export func nextLine(s as Source) {
    return fs.readLine($s.file);
}

/**
 * Close a source and release its child process and temporary files.
 *
 * Returns the decompressor's exit code, so a truncated or corrupt `.gz` can be
 * detected. A stream abandoned early kills the decompressor with `SIGPIPE`, a
 * non-zero code to ignore; use `closeStrict` only on a stream read to the end.
 *
 * @param s {Source} the source to close
 * @return {int} the child's exit code, 0 when there is no child
 * @throws {Error} kind "ngs.source" when the handle could not be closed (closed
 *   twice, say); the child and temporary files are released first
 */
export func close(s as Source) {
    # Release everything before reporting a failed handle close. Letting the
    # `fs.close` error propagate straight out would skip the child reap and the
    # FIFO removal below, leaking both - which is exactly the situation a close
    # is meant to resolve.
    def failure as string init "";
    try {
        fs.close($s.file);
    } catch (e) {
        $failure = $e.message;
    }
    def code as int init 0;
    if ($s.piped) {
        def r as os.Result init os.wait($s.proc);
        $code = $r.exitCode;
        os.release($s.proc);
    }
    discard($s.helper);
    if (len($failure) > 0) {
        core.fail("source", "closing " + $s.origin + " failed (already closed?): " + $failure);
    }
    return $code;
}

/**
 * Close a fully-consumed source, raising if its decompressor failed.
 *
 * @param s {Source} the source to close, already read to exhaustion
 * @throws {Error} kind "ngs.source" when the decompressor exited non-zero
 *   (truncated or corrupt input) or the handle could not be closed
 */
export func closeStrict(s as Source) {
    def code as int init close($s);
    if ($code != 0) {
        core.fail(
            "source",
            "decompressor for " + $s.origin + " exited with code " + convert.toString($code));
    }
}

/**
 * Open a file for streaming line writes, gzip-compressing when the name ends
 * in `.gz`.
 *
 * @param path {string} file to create or overwrite
 * @return {Sink} an open sink; close it with `closeSink`
 * @throws {Error} kind "ngs.source" for a `.gz` target without the `exec`
 *   capability or when `mkfifo` fails; kind "runtime" when the file cannot be
 *   created
 */
export func create(path as string) {
    if (strings.endsWith($path, ".gz")) {
        return createGzip($path);
    }
    return createPlain($path);
}

/**
 * Open an uncompressed file for streaming line writes.
 *
 * @param path {string} file to create or overwrite
 * @return {Sink} an open sink
 * @throws {Error} kind "runtime" when the file cannot be created
 */
export func createPlain(path as string) {
    def f as fs.File init fs.open($path, "write");
    def p as os.Process;
    return Sink{file: $f, proc: $p, helper: "", piped: false, origin: $path};
}

/**
 * Open a gzip-compressed file for streaming line writes at constant memory.
 *
 * Mirrors `openGzip`: lines go into a FIFO that an external `gzip -c`
 * compresses into the target. `closeSink` returns gzip's exit code.
 *
 * @param path {string} `.gz` file to create or overwrite
 * @return {Sink} an open sink
 * @throws {Error} kind "ngs.source" without the `exec` capability, or when
 *   `mkfifo` fails
 */
export func createGzip(path as string) {
    if (not meta.hasCapability("exec")) {
        core.fail(
            "source",
            "streaming gzip output needs the `exec` capability (default jennifer binary)");
    }
    def work as string init makeFifoDir("ngswrite");
    errdefer discard($work);
    def fifo as string init fifoIn($work);
    def cmd as string init "exec gzip -c > " + core.shellQuote($path) + " < " +
        core.shellQuote($fifo);
    def child as os.Process init os.spawn(["sh", "-c", $cmd]);
    errdefer reap($child);
    def f as fs.File init fs.open($fifo, "write");
    return Sink{file: $f, proc: $child, helper: $work, piped: true, origin: $path};
}

/**
 * Write text to a sink verbatim, with no line terminator added.
 *
 * @param k {Sink} an open sink
 * @param text {string} text to write
 */
export func writeText(k as Sink, text as string) {
    fs.writeString($k.file, $text);
}

/**
 * Write one line to a sink, appending a newline.
 *
 * @param k {Sink} an open sink
 * @param line {string} line content without its terminator
 */
export func writeLine(k as Sink, line as string) {
    fs.writeString($k.file, $line + "\n");
}

/**
 * Flush and close a sink, then wait for its compressor to finish.
 *
 * Check the return value: non-zero means the `.gz` output is incomplete.
 *
 * @param k {Sink} the sink to close
 * @return {int} the child's exit code, 0 when there is no child
 * @throws {Error} kind "ngs.source" when the handle could not be closed; the
 *   child and private directory are released first
 */
export func closeSink(k as Sink) {
    # Same ordering as `close`: drain and clean up first, report afterwards.
    def failure as string init "";
    try {
        fs.close($k.file);
    } catch (e) {
        $failure = $e.message;
    }
    def code as int init 0;
    if ($k.piped) {
        def r as os.Result init os.wait($k.proc);
        $code = $r.exitCode;
        os.release($k.proc);
    }
    discard($k.helper);
    if (len($failure) > 0) {
        core.fail("source", "closing " + $k.origin + " failed (already closed?): " + $failure);
    }
    return $code;
}

/**
 * Create a private directory holding a fresh FIFO named `pipe`.
 *
 * The FIFO lives inside a `0700` directory, not directly in the shared temp
 * directory. In `/tmp` it would take reserving a name, deleting it, then
 * calling `mkfifo` on it - a window in which another local user could take the
 * path and hand us a pipe they control. Only the owner can create an entry
 * inside a `0700` directory, so that race does not exist.
 *
 * @param tag {string} directory name prefix
 * @return {string} path of the private directory
 * @throws {Error} kind "ngs.source" when `mkfifo` fails; the `fs` error when the
 *   directory cannot be created
 */
func makeFifoDir(tag as string) {
    def dir as string init fs.makeTempDir("", $tag);
    errdefer discard($dir);
    def made as os.Result init os.run(["mkfifo", fifoIn($dir)]);
    if ($made.exitCode != 0) {
        core.fail("source", "mkfifo failed in " + $dir + ": " + $made.stderr);
    }
    return $dir;
}

/**
 * The FIFO inside a directory made by `makeFifoDir`.
 *
 * @param dir {string} private directory
 * @return {string} path of its FIFO
 */
func fifoIn(dir as string) {
    return path.join($dir, "pipe");
}

/**
 * Remove a temporary file or directory, tolerating one that is already gone.
 *
 * @param target {string} path to remove; "" does nothing
 */
func discard(target as string) {
    if (len($target) == 0 or not fs.exists($target)) {
        return;
    }
    if (fs.isDir($target)) {
        fs.removeAll($target);
        return;
    }
    fs.remove($target);
}

/**
 * Stop a decompressor and drop its handle, on a path that is already failing.
 *
 * Best-effort throughout: it runs from an `errdefer` while another error is
 * propagating, and must not replace it with one of its own.
 *
 * @param child {os.Process} process to stop
 * @return {int} the child's exit code, -1 when it could not be reaped
 */
func reap(child as os.Process) {
    try {
        if (not os.poll($child)) {
            os.kill($child);
        }
        def finished as os.Result init os.wait($child);
        os.release($child);
        return $finished.exitCode;
    } catch (e) {
        return -1;
    }
}
