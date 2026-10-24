# Streaming large files

Every reader in the deck pulls lines from a `source.Source`, and every writer
pushes lines into a `source.Sink`. Both hold one line at a time. The
[format modules](formats.md) wrap `source`; use it directly for line-oriented
files they do not cover.

## Reading lines

```jennifer
use io;
import "@deckyard/ngs/source.j" as source;

def s as source.Source init source.open("big.txt.gz");
while (source.hasNext($s)) {
    io.printf("%s\n", source.nextLine($s));
}
source.closeStrict($s);
```

`source.open` checks the first two bytes for the gzip magic number (`1f 8b`,
RFC 1952) rather than the extension, so a gzipped `.fq` and a plain `.gz` are
both read correctly. `source.openPlain` and `source.openGzip` force one or the
other.

`source.openStdin()` reads standard input through `/dev/stdin` (Unix only):

```sh
samtools view -h aln.bam | jennifer run myscript.j -
```

## Writing lines

```jennifer
def k as source.Sink init source.create("out.bed.gz");
source.writeLine($k, "chr1\t100\t200");
def code as int init source.closeSink($k);
```

`source.create` compresses when the name ends in `.gz`; an output file has no
content to sniff. `source.createPlain` and `source.createGzip` force one or the
other.

## How gzip streaming works

The `compress` library has a streaming compressor (`compress.stream`, `update`,
`finalize`) but no streaming decompressor: `compress.unpack` takes and returns
whole `bytes`. `source.openGzip` therefore:

1. creates a private `0700` directory with `fs.makeTempDir` and runs `mkfifo`
   inside it;
2. spawns `gzip -dc` through `sh -c`, with stdout redirected into the FIFO, so
   nothing accumulates in the interpreter's capture buffer;
3. opens the FIFO with `fs.open` and reads it as an ordinary line handle.

Memory stays at one line regardless of file size, and decompression runs in C.
`source.createGzip` does the same for output, with `gzip -c` reading the FIFO.

The FIFO goes inside its own `0700` directory because creating it directly in
the shared temp directory would mean reserving a name, deleting it and
re-creating it as a FIFO: in that window another local user could take the
path and hand over a pipe they control. Only the owner can create entries in a
`0700` directory.

The path is quoted with `core.shellQuote` before it reaches `sh -c`, so a
filename containing a space, a quote, a backtick, `$(...)` or `;` cannot break
out into shell syntax.

If any step after the FIFO is created fails, an `errdefer` removes the
directory and reaps the decompressor, so a failed open leaves nothing behind,
including when the failure happens inside a caller such as `qc.analyze`.

### Closing, and the exit code

`source.close` closes the read handle (a still-running decompressor then exits
on `SIGPIPE`), waits for the child, removes the FIFO directory and returns the
child's exit code, or `0` when no child was involved. If closing the handle
itself fails, the child and the temporary files are released first and then
`ngs.source` is raised.

A truncated or corrupt `.gz` otherwise looks like an early end of stream.
`source.closeStrict` raises `ngs.source` when the code is non-zero:

```jennifer
source.closeStrict($s);   # raises if the decompressor failed
```

Use `closeStrict` after reading to the end. Use `close` after stopping early:
the `SIGPIPE` exit code is then expected and means nothing went wrong.

### Without the exec capability

`source.openGzipBuffered` decompresses the whole file in memory, writes it to a
temporary file and reads that. It needs no external process, so it runs on
`jennifer-tiny`, but memory is bounded by the decompressed size: use it for
small inputs only.

The deck never falls back to it silently. Without `exec`, `source.openGzip`
raises `ngs.source` naming the fallback, and `source.createGzip` raises
`ngs.source` as well.

## Value semantics and handles

Jennifer copies structs on assignment and argument passing. A `Source` stays
usable after a copy because its `fs.File` and `os.Process` fields are registry
handles: every copy shares the same stream and file position.

The flip side is that a reader cannot carry mutable state such as a pushed-back
line, since a change made inside a called function is not visible to the
caller. That decides the [cursor shapes](formats.md#two-cursor-shapes).
