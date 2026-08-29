# Installing and importing

## Installing

```yaml
# deck.yaml
decks:
  "@jennifer/ngs": "^0.1.0"
```

Two details are worth knowing, because both are easy to get wrong:

- The `@` in `@jennifer/ngs` is *replaced* by the vendor root, so the directory
  is `vendor/jennifer/ngs`, not `vendor/@jennifer/ngs`.
- **Only the deck's `src/` directory is installed**, and its contents land
  directly in `vendor/jennifer/ngs/`. The repository's `docs/`, `examples/`,
  `tools/` and `deck.yaml` are not vendored, so a consumer carries the modules
  and nothing else.

The release archive wraps that `src/` subtree in a single directory, as the
deck spec requires, so vendoring it by hand means stripping those two levels:

```sh
mkdir -p vendor/jennifer/ngs
tar -xzf jennifer-ngs-0.1.0.tar.gz -C vendor/jennifer/ngs \
    --strip-components=2 jennifer-ngs-0.1.0/src
```

`jvc` does this for you: it vendors `src/` and ignores everything beside it.

The vendor root is resolved as `--vendor`, then `JENNIFER_VENDOR`, then the
nearest `vendor/` directory above your program.

## Importing

The entry module carries the deck's identity and its host check. The working
surface lives in the submodules, which you import directly:

```jennifer
import "@jennifer/ngs/" as ngs;            # the entry module
import "@jennifer/ngs/fastq.j" as fastq;
import "@jennifer/ngs/qc.j" as qc;
import "@jennifer/ngs/intervals.j" as intervals;
```

The trailing `/` on the first form expands to the package-named entry module,
`ngs/ngs.j`.

Every module is independent — importing `fastq.j` does not drag in the pipeline
layer. Import what you use.

## Repository layout

```
deck.yaml           the manifest the registry reads
LICENSE             LGPL-3.0-only
src/                the only directory a consumer vendors
  ngs.j             the entry module, named after the deck
  fastq.j ...       the working modules
  *_test.j          white-box test overlays
docs/               this manual
examples/           runnable programs
tools/              CI helpers
```

Every file in `src/` carries an SPDX header and a
`# pragma-jennifer-version: >=0.25.0` floor, which the interpreter checks when
it reads the file — so an interpreter older than that refuses to run the deck
rather than failing somewhere less obvious. `deck.yaml` declares the same floor
under `engines`.

## What your host needs

Two features of the deck shell out, so they need the `exec` capability — the
default `jennifer` binary, not `jennifer-tiny` — and a Unix host with `gzip` and
`mkfifo` on `PATH`:

- streaming gzip input and output ([why](streaming.md#how-gzip-streaming-works))
- [pipeline orchestration](pipelines.md)

Everything else — parsing, interval arithmetic, formatting, QC of uncompressed
input — is pure Jennifer and runs on both binaries.

Check once at startup rather than discovering it partway through a long run:

```jennifer
use io;
import "@jennifer/ngs/" as ngs;

def caps as ngs.Capabilities init ngs.capabilities();
io.printf("%s\n", ngs.describe());
# @jennifer/ngs 0.1.0 (streaming gzip, external tools)

if (not $caps.streamingGzip) {
    io.eprintf("no exec capability: gzip input will have to be buffered\n");
}
```

`ngs.requireFullHost()` is the one-line version that raises instead of
reporting, for a program that cannot do anything useful without both.

## Errors

Every error this deck raises is an `Error` whose `kind` is namespaced
`ngs.<module>` — `ngs.fastq`, `ngs.source`, `ngs.pipeline`, and so on. That
makes them easy to separate from interpreter errors and from each other.

This holds for malformed input too: a non-numeric coordinate raises
`ngs.bed` naming the column, not a bare `runtime` error naming only the text.
Parsing a hostile file is therefore something you can `catch` on kind alone:

```jennifer
try {
    def rec as fastq.Read init fastq.next($r);
} catch (e) {
    io.eprintf("%s: %s\n", $e.kind, $e.message);
}
```

## Limits and hostile input

The deck parses files it did not write, so it is worth being explicit about what
it does and does not defend against.

**What it handles.** Malformed records raise a catchable `ngs.*` error naming the
column, never an uncaught crash. Filenames containing spaces, quotes, backticks,
`$(...)`, `;` or `|` are shell-quoted before reaching `gzip`, so they cannot
inject commands. A failed read removes its FIFO and reaps its decompressor
rather than leaking them. FIFOs are created inside a `0700` directory, so no
other local user can substitute one. `trim` refuses to write its output over its
own input.

**What it does not.** These are inherent to the job, not oversights:

| limit | consequence | what to do |
| --- | --- | --- |
| `source.openGzipBuffered` decompresses whole | a small crafted `.gz` can expand to gigabytes | use the streaming `open`, which is the default |
| a line is read into memory whole | a file with one enormous line costs that much memory | none; line-oriented formats assume sane lines |
| `qc` buffers a batch | bounded by `maxBatchBases`, not by the file | lower it for very long reads |
| `bed.readAll` / `gff.readAll` load everything | proportional to the annotation, by design | use the `next` cursor for large inputs |
| `pipeline` captures stdout and stderr | a tool writing gigabytes to stdout is held in memory | redirect to a file with `pipeline.shell` |
| `pipeline.shell` runs a shell | interpolating untrusted text into it is injection | use `pipeline.job` with an argv, or `core.shellQuote` |

## Running the examples

The repository's `examples/` directory holds runnable programs that import from
`../src/`, so they work straight from a checkout with no vendor tree. They sit
outside `src/` because they are programs rather than modules, and a consumer
should not have to vendor them:

```sh
jennifer run examples/qcreport.j  reads.fastq.gz 200000
jennifer run examples/trimreads.j in.fastq.gz out.fastq.gz
jennifer run examples/intersect.j peaks.bed genes.gff3
jennifer run examples/samfilter.j alignments.sam 30 out.bed
jennifer run examples/vcffilter.j calls.vcf.gz 30
jennifer run examples/pipeline.j  sampleA.fastq.gz sampleB.fastq.gz
```

## Running the tests

Each module has a co-located white-box test overlay in `src/`:

```sh
for suite in src/*_test.j; do jennifer test "$suite"; done
```

A `MODULE_test.j` overlay sees the module's internals, so it calls exported and
private functions alike by their bare names. It also inherits the module's
imports — which is why an overlay must not re-import a module the module under
test already imported, or the alias collides.
