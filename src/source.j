# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Streaming line sources and sinks for large NGS files.
 *
 * The deck's rule is that a 50 GB FASTQ never lands in memory, so every reader
 * in this deck pulls one line at a time from a `Source` and every writer pushes
 * one line at a time into a `Sink`.
 *
 * Plain files stream through an `fs.File` handle. Gzipped files stream through
 * an external `gzip` process whose output is redirected into a FIFO that this
 * module opens as an ordinary read handle - `compress` offers a streaming
 * *compressor* but no streaming *decompressor*, so a pure-Jennifer gzip reader
 * would have to hold the whole decompressed file at once. The FIFO route keeps
 * memory flat at one line regardless of file size, at the cost of needing the
 * `exec` capability (the default `jennifer` binary) and a Unix host.
 *
 * `openGzipBuffered` is the fallback for hosts without `exec`: correct, but it
 * holds the decompressed file in memory, so it is for small inputs only.
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
 * An open line source. Treat it as an opaque handle and pass it around by value:
 * the `fs.File` and `os.Process` fields are registry handles, so copies share
 * the same underlying stream.
 *
 * @field file {fs.File} the readable handle lines are pulled from
 * @field proc {os.Process} the decompressor child, when `piped` is true
 * @field helper {string} private directory or temp file to remove on close
 * @field piped {bool} whether a child process feeds this source
 * @field origin {string} the path this source was opened from, for messages
 */
export def struct Source {
    file as fs.File,
    proc as os.Process,
    helper as string,
    piped as bool,
    origin as string
};

/**
 * An open line sink, the write-side counterpart of `Source`.
 *
 * @field file {fs.File} the writable handle lines are pushed into
 * @field proc {os.Process} the compressor child, when `piped` is true
 * @field helper {string} private directory to remove on close, "" when none
 * @field piped {bool} whether a child process drains this sink
 * @field origin {string} the path this sink was opened for, for messages
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
 * Sniffs the content rather than trusting the extension, so a gzipped `.fq`
 * and a plain `.gz` are both handled correctly.
 *
 * @param path {string} file to inspect
 * @return {bool} true when the first two bytes are 1f 8b
 * @throws {Error} when the file cannot be opened
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
 * @return {Source} an open source; close it with `close`
 * @throws {Error} when the file is missing, or gzipped without `exec` available
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
 * @throws {Error} when the file is missing
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
 * Open standard input as a line source, so the deck composes with Unix pipes.
 *
 * Uses `/dev/stdin`, so it is Unix-only.
 *
 * @return {Source} an open source reading the process's stdin
 * @throws {Error} when /dev/stdin is not available
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
 * Spawns `gzip -dc` with its stdout redirected into a private FIFO, then reads
 * the FIFO as a normal line handle. Memory stays flat no matter how large the
 * file is.
 *
 * @param path {string} gzipped file to read
 * @return {Source} an open source
 * @throws {Error} without the `exec` capability, or when `mkfifo` fails
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
 * Holds the entire decompressed file in memory and then on disk as a temp file,
 * so use it only for small inputs; `openGzip` is the streaming path.
 *
 * @param path {string} gzipped file to read
 * @return {Source} an open source over a temporary decompressed copy
 * @throws {Error} when the file is missing or not valid gzip
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
 * Whether another line can be read from the source.
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
 * @throws {Error} when the source is already exhausted; guard with `hasNext`
 */
export func nextLine(s as Source) {
    return fs.readLine($s.file);
}

/**
 * Close a source and release its child process and temporary files.
 *
 * Returns the decompressor's exit code so a truncated or corrupt `.gz` can be
 * detected. Note that abandoning a stream early makes the decompressor die of
 * `SIGPIPE`, which is a non-zero code you should ignore; use `closeStrict` only
 * when the stream was consumed to the end.
 *
 * @param s {Source} the source to close
 * @return {int} the child's exit code, or 0 when no child was involved
 * @throws {Error} kind `ngs.source` when the handle could not be closed - after
 *   the child and temporary files have been released regardless
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
 * @throws {Error} when the decompressor exited non-zero (truncated/corrupt input)
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
 * @throws {Error} for a `.gz` target without the `exec` capability
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
 */
export func createPlain(path as string) {
    def f as fs.File init fs.open($path, "write");
    def p as os.Process;
    return Sink{file: $f, proc: $p, helper: "", piped: false, origin: $path};
}

/**
 * Open a gzip-compressed file for streaming line writes at constant memory.
 *
 * Mirrors `openGzip`: lines are written into a FIFO that an external `gzip`
 * process compresses into the target file.
 *
 * @param path {string} `.gz` file to create or overwrite
 * @return {Sink} an open sink
 * @throws {Error} without the `exec` capability, or when `mkfifo` fails
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
 * Write text to a sink verbatim, adding no line terminator.
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
 * Close a sink, flush it and wait for its compressor to finish.
 *
 * @param k {Sink} the sink to close
 * @return {int} the child's exit code, or 0 when no child was involved
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
 * Create a private directory holding a fresh FIFO, and return the directory.
 *
 * The FIFO lives inside a `0700` directory rather than directly in the shared
 * temp directory. Putting it in `/tmp` would mean reserving a name, deleting it,
 * then calling `mkfifo` on it - a window in which another local user could take
 * the path and hand us a pipe they control. Only the owner can create an entry
 * inside a `0700` directory, so that race does not exist.
 *
 * @param tag {string} prefix used in the directory name, for debuggability
 * @return {string} path to a private directory containing `pipe`
 * @throws {Error} when the directory or the FIFO cannot be created
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
 * @param dir {string} the private directory
 * @return {string} path to its FIFO
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
 * Every step is best-effort: this runs from an `errdefer` while another error is
 * propagating, and must not replace it with one of its own.
 *
 * @param child {os.Process} the process to stop
 * @return {int} the child's exit code, or -1 when it could not be reaped
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
