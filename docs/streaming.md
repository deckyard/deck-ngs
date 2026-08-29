# Streaming large files

The deck's central constraint is that a 50 GB FASTQ never lands in memory. Every
reader here is a pull cursor over a `source.Source`, and every writer pushes into
a `source.Sink`. Both hold one line at a time.

You rarely use `source` directly — [the format modules](formats.md) wrap it —
but understanding it explains the shape of everything above it.

## Reading lines

```jennifer
use io;
import "@jennifer/ngs/source.j" as source;

def s as source.Source init source.open("big.txt.gz");
while (source.hasNext($s)) {
    io.printf("%s\n", source.nextLine($s));
}
source.closeStrict($s);
```

`source.open` sniffs the file's first two bytes for the gzip magic number rather
than trusting the extension, so a gzipped `.fq` and a plain `.gz` are both
handled correctly. `source.openPlain` and `source.openGzip` force one or the
other.

`source.openStdin()` reads the process's standard input via `/dev/stdin`, which
is how the deck composes with Unix pipes:

```sh
samtools view -h aln.bam | jennifer run myscript.j -
```

## Writing lines

```jennifer
def k as source.Sink init source.create("out.bed.gz");
source.writeLine($k, "chr1\t100\t200");
def code as int init source.closeSink($k);
```

`source.create` compresses when the target name ends in `.gz` — by extension
here, since an output file has no content to sniff yet.

## How gzip streaming works

The `compress` library ships a streaming *compressor* — `compress.stream`,
`update`, `finalize` — but no streaming *decompressor*. `compress.unpack` takes
and returns whole `bytes`, so decompressing a 50 GB FASTQ through it would need
50 GB of memory, which defeats the purpose.

So `source.openGzip` does this instead:

1. create a private `0700` directory with `fs.makeTempDir` and `mkfifo` a FIFO
   *inside* it — putting the FIFO straight in the shared temp directory would
   mean reserving a name, deleting it and then creating it, a window in which
   another local user could take the path and hand us a pipe they control;
2. `os.spawn` a `gzip -dc` whose stdout the shell redirects into that FIFO —
   redirected, so nothing accumulates in the interpreter's capture buffer;
3. `fs.open` the FIFO and read it as an ordinary line handle.

Memory stays flat at one line no matter how large the file is, and the
decompression itself runs in C. `source.createGzip` mirrors it for output, with
the FIFO on the writing side.

The path is quoted with `core.shellQuote` before it reaches `sh -c`, so a
filename containing a space, a quote, a backtick, a `$(...)` or a semicolon
cannot break out into shell syntax. The deck's test suite reads files with all
of those in their names and checks that nothing else executes.

If any step after the FIFO is created fails, an `errdefer` removes the directory
and reaps the decompressor, so a failed open leaves nothing behind — including
when the failure happens deep in a caller like `qc.analyze`.

### Closing, and why it returns a code

`source.close` closes the read handle first — which makes the decompressor exit
on `SIGPIPE` if it is still running — then waits for the child, removes the
FIFO, and returns the child's exit code.

That code is worth checking, because a **truncated or corrupt `.gz` is otherwise
silent**: the reader simply sees end-of-stream early. `source.closeStrict` is the
checking form:

```jennifer
source.closeStrict($s);   # raises if the decompressor failed
```

Use `closeStrict` when you consumed the stream to the end. Use plain `close`
when you stopped early — abandoning a stream gives the decompressor a `SIGPIPE`
and a non-zero code that means nothing went wrong.

### Without the exec capability

`source.openGzipBuffered` decompresses the whole file into memory and then to a
temporary file, and reads that. It needs no external process, so it works on
`jennifer-tiny` — but it is bounded by the decompressed size, so it is for small
inputs only.

The deck never falls back to it silently. `source.openGzip` raises a clear error
naming the fallback, because quietly turning a constant-memory read into a 50 GB
allocation is worse than failing.

## Value semantics and handles

A `Source` is a struct, and Jennifer copies structs on assignment and argument
passing. That is safe here because the fields that matter — the `fs.File` and
the `os.Process` — are registry handles: copies share the same underlying
stream. Passing a `Source` into a function does not fork the file position.

The consequence to know is the other direction: a reader **cannot carry mutable
state**, because a mutation inside a called function would not be visible to the
caller. That is why the record cursors are shaped the way they are — see
[the note on end-of-stream](formats.md#two-cursor-shapes).
