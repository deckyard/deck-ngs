# Installing and importing

## Installing

```yaml
# deck.yaml
decks:
  "@deckyard/ngs": "^0.1.0"
```

The `@` in `@deckyard/ngs` is replaced by the vendor root, so the deck lands in
`vendor/deckyard/ngs/`, not `vendor/@deckyard/ngs/`. Only the deck's `src/` is
vendored, and its contents go directly into that directory; `docs/`,
`examples/`, `tools/` and `deck.yaml` stay behind.

The release archive, `deckyard-ngs-<version>.tar.gz`, wraps `src/` in one
top-level directory. To vendor it by hand, strip both levels:

```sh
mkdir -p vendor/deckyard/ngs
tar -xzf deckyard-ngs-0.1.0.tar.gz -C vendor/deckyard/ngs \
    --strip-components=2 deckyard-ngs-0.1.0/src
```

The vendor root is resolved as `--vendor`, then `JENNIFER_VENDOR`, then the
nearest `vendor/` directory above the program.

Every file in `src/` carries a `# pragma-jennifer-version: >=0.25.0` header, so
an older interpreter refuses to load the deck. `deck.yaml` declares the same
floor under `engines`.

## Importing

```jennifer
import "@deckyard/ngs/" as ngs;            # the entry module, ngs/ngs.j
import "@deckyard/ngs/fastq.j" as fastq;
import "@deckyard/ngs/qc.j" as qc;
import "@deckyard/ngs/intervals.j" as intervals;
```

The trailing `/` expands to the entry module, which holds the deck's name,
version and host check. The working surface is in the submodules; import the
ones you use.

## What the host needs

Two features run external processes. They need the `exec` capability (the
default `jennifer` binary, not `jennifer-tiny`) and a Unix host with `gzip` and
`mkfifo` on `PATH`:

- streaming gzip input and output ([how it works](streaming.md#how-gzip-streaming-works))
- [pipeline orchestration](pipelines.md)

Parsing, interval arithmetic, formatting and QC of uncompressed input are pure
Jennifer and run on both binaries.

Check at startup:

```jennifer
use io;
import "@deckyard/ngs/" as ngs;

def caps as ngs.Capabilities init ngs.capabilities();
io.printf("%s\n", ngs.describe());
# @deckyard/ngs 0.1.0 (streaming gzip, external tools)

if (not $caps.streamingGzip) {
    io.eprintf("no exec capability: gzip input will have to be buffered\n");
}
```

`ngs.requireFullHost()` raises `ngs.ngs` instead when `exec` is missing.

## Errors

Every error the deck raises has `kind` set to `ngs.<module>`: `ngs.fastq`,
`ngs.source`, `ngs.pipeline` and so on. Malformed input is included: a
non-numeric coordinate raises `ngs.bed` naming the record, not a bare `runtime`
error. A reader can catch on kind alone:

```jennifer
try {
    def rec as fastq.Read init fastq.next($r);
} catch (e) {
    io.eprintf("%s: %s\n", $e.kind, $e.message);
}
```

The exception is a failure inside the interpreter's own file calls, such as
`source.openStdin` on a standard input that cannot be reopened (a socket):
that surfaces as `runtime`.

## Limits and hostile input

The deck defends against:

- malformed records, which raise a catchable `ngs.*` error naming the record;
- filenames containing spaces, quotes, backticks, `$(...)`, `;` or `|`, which
  are shell-quoted before they reach `gzip`;
- a failed open, which removes its FIFO and reaps its decompressor;
- another local user substituting the FIFO, which is created inside a `0700`
  directory;
- `trim` writing its output over its own input, which it refuses.

It does not defend against these, which come with the job:

| limit | consequence | what to do |
| --- | --- | --- |
| `source.openGzipBuffered` decompresses whole | a small crafted `.gz` can expand to gigabytes | use the streaming `source.open`, the default |
| a line is read into memory whole | a file with one enormous line costs that much memory | none; line formats assume sane lines |
| `qc` buffers a batch | bounded by `maxBatchBases`, not by the file | lower it for long reads |
| `bed.readAll` / `gff.readAll` load everything | memory proportional to the annotation | use the `next` cursor for large inputs |
| `pipeline` captures stdout and stderr | a tool writing gigabytes to stdout is held in memory | redirect to a file with `pipeline.shell` |
| `pipeline.shell` runs a shell | interpolating untrusted text into it is command injection | use `pipeline.job` with an argv, or `core.shellQuote` |

## Examples

`examples/` in the repository holds runnable programs. They import from
`../src/`, so they run from a checkout without a vendor tree:

```sh
jennifer run examples/qcreport.j  reads.fastq.gz 200000
jennifer run examples/trimreads.j in.fastq.gz out.fastq.gz
jennifer run examples/intersect.j peaks.bed genes.gff3
jennifer run examples/samfilter.j alignments.sam 30 out.bed
jennifer run examples/vcffilter.j calls.vcf.gz 30
jennifer run examples/pipeline.j  sampleA.fastq.gz sampleB.fastq.gz
```

Running the tests and the other checks is covered in
[Working on the deck](development.md).
