# Design

## Scope

The deck is the glue and light-I/O layer for next-generation sequencing
data:

- **Streaming format I/O** for FASTQ (plain or gzipped) and the
  tab-delimited SAM, VCF, BED and GFF3/GTF formats: read, filter, convert
  and write record by record, so a 50 GB file never lands in memory.
- **Read QC and preprocessing** (FastQC-lite, fastp-lite): per-base quality,
  length and GC distributions, adapter detection, quality and adapter
  trimming, length filtering and subsampling, each in one streaming pass.
- **Interval arithmetic** (bedtools-lite): sort, merge, intersect, subtract,
  coverage and complement over BED and GFF features.
- **Pipeline orchestration** (Snakemake-lite): run external tools as jobs,
  skip jobs whose outputs are up to date, fan out over samples in parallel,
  manage intermediate files.

NGS breaks the model of the pure-Jennifer decks on two axes. Files are huge
(FASTQ.gz in the tens of GB, BAM in the hundreds), so everything streams.
And the heavy compute (alignment, assembly, variant calling) is
SIMD-optimised C that an interpreter is 100 to 1000 times too slow to
replace. So the deck deliberately stays out of:

- read **alignment** (BWA, Bowtie2, minimap2), de-novo **assembly**
  (SPAdes) and **variant calling** (GATK): drive the native tools through
  `pipeline`;
- **BAM and CRAM**: they also need BGZF (blocked gzip with random access),
  which `compress` does not provide, and have little standalone value
  without the downstream compute the deck omits.

## Decisions

- **A modules deck.** Consumers import the submodule they need
  (`import "@deckyard/ngs/qc.j" as qc;`). The entry module carries only
  identity (`VERSION`, `NAME`) and the host probe (`capabilities`,
  `describe`, `requireFullHost`). Every module imports `core`; the readers
  import `source`; the format modules import `intervals`.
- **Pure Jennifer, with gzip through a subprocess.** The parsers are
  Jennifer code over `fs` line reads. `compress` has no streaming
  decompressor, so a gzipped input is read by running `gzip -dc` into a FIFO
  and reading the FIFO as an ordinary file; gzipped output goes the other
  way through `gzip -c`. `openGzipBuffered` is the fallback for hosts
  without `exec`: it decompresses the whole file in memory and is only for
  small inputs.
- **`exec` is checked at run time.** The manifest declares the `exec`
  capability, but the modules test `meta.hasCapability("exec")` where they
  need it instead of declaring a file pragma, so parsing, intervals and QC
  still load on a host without it.
- **Gzip by content when reading, by name when writing.** `source.open`
  sniffs the gzip magic bytes, so a gzipped file named `.fq` reads fine;
  `source.create` picks gzip by the `.gz` suffix.
- **Two reader shapes.** `fastq.Reader` uses `hasNext` and `next` over
  fixed four-line records. The BED, GFF, SAM and VCF readers return a `Next`
  enum (a record or `Done`), because skipping comment and header lines can
  consume the last line, which a `hasNext` peek cannot express.
- **Coordinates.** `intervals.Interval` is 0-based half-open, like BED.
  GFF, SAM and VCF are 1-based; their `toInterval` converts. Intervals are
  built with `make` or `named`, which validate; strands are `+`, `-` or `.`,
  so a GFF3 `?` strand is an error when converted.
- **Interval order.** Sort keys are zero-padded strings of `KEY_WIDTH`
  digits, and a coordinate at or past `MAX_COORDINATE` (10^12) is refused at
  construction so it cannot silently missort. Chromosome order is
  lexicographic (`chr10` before `chr2`), matching `sort -k1,1 -k2,2n`.
  Already sorted input is returned unchanged after one linear check.
- **Linear interval sweeps.** Intersect and overlap counting keep an active
  list of targets spanning the query start, pruned per query, so one long
  target cannot make the sweep quadratic. Coverage and subtraction sweep
  merged, disjoint targets with a single cursor.
- **Strict records.** Parsers validate each record (field counts, numbers,
  the quality alphabet) and report the record that failed. Writers refuse
  tab, CR and LF inside fields. Commas are never escaped in GFF3 attributes
  or VCF INFO values, because they are the multi-value separator; repeated
  GFF/GTF keys are joined with `,`, repeated VCF INFO keys are an error.
- **Adapter search.** An exact substring match first, then a seed search on
  the first bases of the adapter for 3' partial adapters. Exact and
  case-sensitive, no mismatches. An empty adapter is an error.
- **Safe outputs.** `trim.run` and the subsampling functions refuse to write
  over their own input, write to a temporary file beside the output and
  rename it over the output on success (removing it on failure). The rename
  also protects a hard-linked output, which a same-file check cannot see.
- **Close semantics.** `source.close` returns the decompressor's exit code
  (a stream abandoned early gets SIGPIPE, which callers ignore);
  `closeStrict` raises on a non-zero code and is for fully read streams. A
  non-zero code when closing a gzip sink means truncated output, and
  `trim.run` raises on it. `qc.analyze` uses the lenient close, so a
  truncated gzip input is not reported there.
- **Pipelines without a shell.** Jobs run through `os.run` with an argv;
  `shell()` is the explicit `sh -c` opt-in. `pipeline.needsRun` treats an
  input whose mtime equals the oldest output as stale, because coarse
  timestamps make ties common and skipping a stale job is worse than
  re-running a current one. `runParallel` runs waves of at most `width`
  jobs; it does no dependency resolution.
- **Performance shapes.** QC counts characters with one native scan per
  batch (`DEFAULT_BATCH` records or `DEFAULT_BATCH_BASES` bases, whichever
  comes first); per-position quality cannot be counted that way and is
  opt-in. Hot readers call `fs` directly rather than through the `source`
  wrappers, per-record loops accumulate into locals rather than a struct,
  and `runParallel` spawns each job from a one-job helper so a spawn does
  not copy the whole job list. Numbers are in the Performance chapter.

## Errors

`core.fail` throws `Error{kind: "ngs.<what>"}`: `ngs.source`, `ngs.fastq`,
`ngs.bed`, `ngs.gff`, `ngs.sam`, `ngs.vcf`, `ngs.intervals`, `ngs.qc`,
`ngs.trim` and `ngs.pipeline`. Record parsers run each record in one `try`
and pass the error to `core.reraise`, which rethrows `ngs.*` errors
unchanged and wraps anything else (such as a `runtime` error from a number
conversion) as `ngs.<what>` with the record appended. Record parsing works
on one line in memory, so no transport, timeout or cap error passes through
that wrapper.

Some paths still surface the interpreter's own `runtime` kind: reading past
the end of a source, a failure to create a plain output file, and opening
stdin when it is a socket. Their docblocks say so. File-system errors
elsewhere are passed on unchanged. The deck talks to no server and has no
error struct.

## Security and limits

- FIFOs live in a private `fs.makeTempDir` directory, which closes the race
  between reserving a name in `/tmp` and creating the FIFO.
- Paths reach the gzip shell command through single-quote escaping, after
  `--`, so a file name cannot inject options or commands. Pipeline jobs take
  an argv and run without a shell unless the caller asks for one.
- `source.open` refuses a non-regular path (a FIFO or piped stdin): sniffing
  reads it and then reopens it, which would drop bytes or block. The plain
  and stdin openers open once and accept them.
- `pipeline.cleanUp` refuses `""`, `.`, `..` and `/` by name, then compares
  real paths and refuses the root, the current directory and any ancestor of
  it.
- `openGzipBuffered` holds the whole decompressed file in memory, so a small
  compressed input can expand without bound. Use it only on trusted, small
  files.
- BED scores beyond 10^15 and coordinates at or past `MAX_COORDINATE` are
  refused.

## Testing and validation

- Every module has an overlay; together they hold 217 tests. Fixtures are
  written at test time into a git-ignored directory in the checkout (or the
  system temporary directory for older tests and for FIFOs, which some
  shared folders refuse). There is no `testdata/`.
- The tests need `sh`, `gzip`, `mkfifo` and `touch` on the host: gzip FIFO
  round trips, pipeline jobs, and exact mtimes staged with `touch -d`.
- `tools/check-examples.sh` builds a tiny FASTQ, BED, GFF3, SAM and VCF set
  and runs five examples on it; `examples/pipeline.j` is only parse-checked,
  because CI runners have no aligners.

Not covered: real aligners and samtools are never run in CI, and the
formats are tested against small hand-written files rather than published
datasets.

## Known limitations and open questions

- `sam.open` reads the header in a separate pass, so a gzipped SAM is
  decompressed twice; `sam.openStdin` leaves the header empty.
- `capabilities` reports `streamingGzip` from the capability alone; a caller
  must also check `gzipOnPath`.
- `ngs.VERSION` is kept by hand and must be bumped with the manifest;
  `tools/check-extra.sh` fails while the two differ.
- `jennifer-tiny` is not a supported engine: it lacks `exec` and the
  capability probe the gzip and pipeline layers need.

## Dependencies on the interpreter

- `strings.count`, the native scan QC counts with; it is why the floor is
  `>=0.26.0`.
- `os.run`, `spawn` and `task.waitAll` for pipelines, `fs` line reads and
  `meta.hasCapability`.
- A `spawn` copies every local in scope; `jennifer`'s planned M30.7
  (capture only the locals a `spawn` body reads) would remove the need for
  the one-job helper in `runParallel`.
