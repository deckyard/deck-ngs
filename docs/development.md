# Working on the deck

This page is for people changing the deck rather than using it. Everything here
runs against a checkout; none of it is needed to consume `@jennifer/ngs`.

## The gates

Six checks stand between a change and a green build. All of them run locally,
and CI runs exactly these, through exactly these scripts:

```sh
jennifer fmt -l src/*.j examples/*.j     # formatting
jennifer lint src/*.j examples/*.j       # lint
tools/check-headers.sh                   # SPDX pair and engine floor, against deck.yaml
tools/check-overlays.sh                  # every module has a co-located overlay
tools/check-cheatsheet.sh                # every export appears in docs/cheatsheet.md
tools/check-examples.sh                  # every example still runs
for suite in src/*_test.j; do jennifer test "$suite"; done
```

`fmt -l` lists unformatted files and exits 1 when there are any; `lint` exits 1
on a warning and 0 when clean. The scripts guard against drift that neither tool
can see:

- **`check-headers.sh`** reads the licence and the engine floor out of
  `deck.yaml` and checks that every `.j` file's header still states both. Those
  values are deliberately declared twice, and this is what keeps the copies
  honest. `jennifer lint` rejects a *malformed* pragma but says nothing about a
  missing one, so without this a new module could land bare and lint clean.
- **`check-overlays.sh`** requires a co-located `MODULE_test.j` beside every
  module, and a module beside every overlay. That is
  [jvc's publish gate](https://github.com/jennifer-language/jvc/blob/main/docs/manifest.md),
  enforced here so it fails at push time rather than at publish time.
- **`check-cheatsheet.sh`** greps every `export` out of `src/` and fails if a
  name is missing from [the cheatsheet](cheatsheet.md). Adding an export without
  documenting it fails CI.
- **`check-examples.sh`** builds a small FASTQ/BED/GFF3/SAM/VCF fixture set in a
  temporary directory and runs each example against it, because every example
  here takes input files. `examples/pipeline.j` drives external aligners a
  runner does not have, so it is parse-checked instead.

These scripts are shared verbatim with the other deck repositories, and
`tools/deck-lib.sh` is what makes that possible: it discovers the manifest
(`deck.toml`, `deck.yaml`, `deck.yml` or `deck.json`), the deck's identity, and
whether it is a flat set of modules or a module plus `include`d topics. Nothing
about this deck is written into a workflow.

The suite is **143 tests across twelve overlays**, covering gzip round-trips
through the FIFO path, interval algebra against known bedtools semantics,
coordinate conversion in both directions, parser edge cases (minimal records,
symbolic VCF alleles, `##FASTA` sections, partial adapters), QC values checked
against an independent reference implementation, and real concurrency in
`runParallel`.

## Building the book

The manual is the Markdown in `docs/`, with `docs/SUMMARY.md` as its outline and
`grimoire.toml` as its configuration. It is built by
[Grimoire](https://grimoire.jennifer-lang.dev/), Jennifer's own static site
generator, which ships as a container - so building the book needs nothing
installed but Docker:

```sh
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" \
    ghcr.io/jennifer-language/grimoire:latest build      # -> site/
```

From a Grimoire checkout instead, with a live preview on `:8080`:

```sh
/path/to/grimoire/bin/grimoire build
/path/to/grimoire/bin/grimoire serve --watch
```

A chapter named in `SUMMARY.md` with no file behind it fails the build, which is
what stops a dead nav entry reaching the published site. The book is published as
a site only - `[pdf] enabled = false` in `grimoire.toml`.

## What CI does

Three workflows, byte-identical to the ones in the other deck repositories:

| workflow | when | what |
| --- | --- | --- |
| [`test.yml`](https://github.com/jennifer-language/deck-ngs/blob/main/.github/workflows/test.yml) | push, PR | fmt, lint, the four guards, every suite, every example |
| [`docs.yml`](https://github.com/jennifer-language/deck-ngs/blob/main/.github/workflows/docs.yml) | push, PR | builds the book; publishes to Pages from `main` |
| [`release.yml`](https://github.com/jennifer-language/deck-ngs/blob/main/.github/workflows/release.yml) | version tag | re-runs the gate, then cuts the release |

Both the tests and the docs run in containers, so a runner needs nothing
installed: `ghcr.io/jennifer-language/jennifer:dev` for the interpreter and
`ghcr.io/jennifer-language/grimoire:latest` for the book. The deck's `>=0.25.0`
floor is not a released version yet and a dev build bypasses it; once 0.25.0
ships, `JENNIFER_IMAGE` becomes `:0.25` and nothing else moves.

`test.yml` also carries a `tiny` job that only does anything when the manifest
lists `jennifer-tiny` under `engines`. This deck does not - streaming gzip and
pipeline orchestration need `exec` - so the job resolves to a no-op here.

Publishing the book needs GitHub Pages set to **Build and deployment → Source →
GitHub Actions** in the repository settings.

## Cutting a release

Tag the commit with the version already in `deck.yaml`. The `release` job refuses
a tag that disagrees with the manifest, and accepts either `0.2.0` or `v0.2.0`.

`tools/pack.sh` builds the archive, and it is the shape
[the deck spec](https://github.com/jennifer-language/jvc/blob/main/docs/deck-spec.md)
describes: `jennifer-ngs-<version>.tar.gz`, one wrapping directory, and a `src/`
subtree inside it, with the manifest, licence and README beside it. `jvc` vendors
`src/` and ignores everything else. Alongside it go two forms of the same digest
- a `sha256sum -c` file for humans, and the `sha256:<hex>` the registry records
as its integrity pin.

Installing that archive by hand means stripping the two leading levels:

```sh
mkdir -p vendor/jennifer/ngs
tar -xzf jennifer-ngs-0.2.0.tar.gz -C vendor/jennifer/ngs \
    --strip-components=2 jennifer-ngs-0.2.0/src
```
