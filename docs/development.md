# Working on the deck

For changing the deck; none of this is needed to use it.

## The gate

`tools/check-all.sh` runs every check CI runs, in CI's order, continues past a
failure and ends with "the gate is green" or the list of failed steps. The
checks one by one:

```sh
tools/check-shared.sh                    # shared files match the template
jennifer fmt -l src/*.j examples/*.j     # formatting
jennifer lint src/*.j examples/*.j       # lint
tools/check-headers.sh                   # SPDX pair and engine floor, against deck.yaml
tools/check-overlays.sh                  # every module has a co-located overlay
tools/check-src.sh                       # src/ holds only .j files
tools/check-coverage.sh                  # every export is called from an overlay
tools/check-docs.sh                      # docs, signatures and docblocks match the code
tools/check-identity.sh                  # no template placeholder left
for suite in src/*_test.j; do jennifer test "$suite"; done
tools/check-examples.sh                  # every example still runs
```

`fmt -l` lists unformatted files and exits 1 if there are any; `lint` exits 1
on a warning and 0 on info. The scripts catch drift neither tool sees:

- `check-headers.sh` reads the licence and the engine floor from `deck.yaml`
  and checks every `.j` header states both. `jennifer lint` rejects a malformed
  pragma but not a missing one.
- `check-overlays.sh` requires a `MODULE_test.j` beside every module and a
  module beside every overlay. This is
  [jvc's publish gate](https://github.com/jennifer-language/jvc/blob/main/docs/manifest.md),
  checked at push time instead of publish time.
- `check-docs.sh` compares [the cheatsheet](cheatsheet.md) and every
  deck-module call in the docs with the code: names, parameter names, arity,
  struct fields, enum variants, constant values and docblock tags.
  `docs/.cheatsheet-todo`, when present, lists tolerated older gaps; it only
  shrinks.
- `check-coverage.sh` fails when an export is never called from an overlay;
  `tools/.coverage-todo` is its backlog. The remaining entry,
  `source.createGzip`, can only be exercised by spawning `gzip`.
- `check-examples.sh` builds a FASTQ/BED/GFF3/SAM/VCF fixture set in
  `.ci-fixtures/` inside the checkout, runs each example against it with
  relative paths, and removes it. `examples/pipeline.j` drives external
  aligners and is only parse-checked.

The `tools/` scripts other than `check-examples.sh` are shared byte-for-byte
with the other deck repositories. `tools/deck-lib.sh` discovers the manifest
(`deck.toml`, `deck.yaml`, `deck.yml` or `deck.json`), the deck's identity and
whether it is a set of modules or a module with `include`d topics.

## Tests

```sh
jennifer test src/qc_test.j              # one suite; there is no finer filter
for suite in src/*_test.j; do jennifer test "$suite"; done
```

There are 172 tests in twelve overlays. They cover gzip round-trips through the
FIFO, interval operations against bedtools semantics, coordinate conversion in
both directions, parser edge cases (minimal records, symbolic VCF alleles,
`##FASTA` sections, partial adapters), QC values against hand-computed counts,
and concurrency in `runParallel`.

`MODULE_test.j` is a white-box overlay of `MODULE.j`: it calls exported and
private functions by their bare names and inherits the module's imports, so
re-importing one of them is an alias collision.

In CI `jennifer` runs in a container that mounts only the checkout. A test that
needs files creates them under `.ci-fixtures/` (git-ignored) and passes
relative paths. Tests that read standard input only open and close it; on a
host where it cannot be reopened (a socket) they check nothing.

## Building the book

The manual is the Markdown in `docs/`, with `docs/SUMMARY.md` as its outline
and `grimoire.toml` as its configuration. It is built by
[Grimoire](https://grimoire.jennifer-lang.dev/), which ships as a container:

```sh
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" \
    ghcr.io/jennifer-language/grimoire:latest build      # -> site/
```

From a Grimoire checkout, with a live preview on `:8080`:

```sh
/path/to/grimoire/bin/grimoire build
/path/to/grimoire/bin/grimoire serve --watch
```

A chapter in `SUMMARY.md` with no file behind it fails the build. The book is
published as a site only (`[pdf] enabled = false` in `grimoire.toml`).

## CI

Three workflows, byte-identical across the deck repositories:

| workflow | when | what |
| --- | --- | --- |
| [`test.yml`](https://github.com/deckyard/deck-ngs/blob/main/.github/workflows/test.yml) | push, PR | the gate above |
| [`docs.yml`](https://github.com/deckyard/deck-ngs/blob/main/.github/workflows/docs.yml) | push, PR | builds the book; publishes to Pages from `main` |
| [`release.yml`](https://github.com/deckyard/deck-ngs/blob/main/.github/workflows/release.yml) | version tag | re-runs the gate, then cuts the release |

Tests and docs run in containers, so a runner needs nothing installed.
`tools/deck-lib.sh` picks the interpreter image from the floor: a floor at or
below the newest release (0.25) runs on that release's image, here
`ghcr.io/jennifer-language/jennifer:0.25`; a newer floor runs on `:dev`, which
bypasses the floor check. `JENNIFER_IMAGE` overrides the choice. The book uses
`ghcr.io/jennifer-language/grimoire:latest`.

`test.yml` has a `tiny` job that runs only when the manifest lists
`jennifer-tiny` under `engines`. This deck does not, since streaming gzip and
`pipeline` need `exec`.

Publishing the book needs GitHub Pages set to Build and deployment > Source >
GitHub Actions in the repository settings.

## Cutting a release

Tag the commit with the version in `deck.yaml`, written plain (`0.2.0`, no `v`).
`release.yml` refuses any other tag.

`tools/pack.sh` builds the archive in the shape
[the deck spec](https://github.com/jennifer-language/jvc/blob/main/docs/deck-spec.md)
describes: `deckyard-ngs-<version>.tar.gz`, one wrapping directory holding
`src/` and, beside it, the manifest, licence and README. `jvc` vendors `src/`
and ignores the rest. Two forms of the digest go with it: a `sha256sum -c` file
and the `sha256:<hex>` the registry records as its integrity pin.

[Installing and importing](getting-started.md#installing) shows how to vendor
the archive by hand.
