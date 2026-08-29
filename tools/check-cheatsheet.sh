#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
#
# Guard against docs/cheatsheet.md drifting from the code.
#
# The cheatsheet is hand-committed rather than generated at build time, so this
# checks the one property that matters: every exported name appears in it.
# Adding an export without documenting it fails CI.
#
# How a name is written follows the deck's shape. A modules deck has one
# namespace per module and qualifies every name with it (`fastq.open(`). A
# topics deck presents a single namespace that the reader aliases as it likes,
# so its cheatsheet writes the names bare (`translate(`).
#
# A deck with no cheatsheet skips the check rather than failing it.
#
# Known gaps go in docs/.cheatsheet-todo, one name per line. Those are reported
# as a backlog instead of failing, so a deck can adopt this guard without a red
# build while still failing the moment a *new* export goes undocumented. Shrink
# the file; do not grow it.
#
# Usage: tools/check-cheatsheet.sh
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck source=tools/deck-lib.sh
source tools/deck-lib.sh
deck_load

ref="docs/cheatsheet.md"
if [[ ! -f "$ref" ]]; then
    echo "no ${ref} in this deck - skipping the export coverage check"
    exit 0
fi

todo="docs/.cheatsheet-todo"
known() {
    [[ -f "$todo" ]] || return 1
    grep -qxF -- "$1" "$todo"
}

missing=0
backlog=0
checked=0
while read -r file; do
    case "$file" in *_test.j) continue ;; esac
    if [[ "$DECK_SHAPE" == topics ]]; then
        prefix=""
    else
        prefix="$(basename "$file" .j)."
    fi

    while read -r kind name; do
        checked=$((checked + 1))
        qualified="${prefix}${name}"
        # A function is cited with its call parenthesis; a type or constant is
        # cited by name, which may sit inside a larger type expression such as
        # `list of Bridge`.
        if [[ "$kind" == "func" ]]; then
            token="${qualified}("
        else
            token="$qualified"
        fi

        if grep -qF -- "$token" "$ref"; then
            continue
        fi
        if known "$qualified"; then
            backlog=$((backlog + 1))
            continue
        fi
        echo "undocumented export: ${qualified} (${kind})" >&2
        missing=$((missing + 1))
    done < <(
        grep -oE '^export (func|def const|def struct|def enum) [A-Za-z][A-Za-z0-9_]*' "$file" \
            | sed -E 's/^export (func|def const|def struct|def enum) /\1 /; s/^def //'
    )
done < <(find src -name '*.j' | sort)

if (( missing > 0 )); then
    echo "" >&2
    echo "$missing of $checked exported names are missing from $ref." >&2
    echo "Document them there, then re-run tools/check-cheatsheet.sh." >&2
    exit 1
fi

if (( backlog > 0 )); then
    echo "all $checked exported names are documented in $ref,"
    echo "  except $backlog carried in ${todo} - shrink that file"
else
    echo "all $checked exported names are documented in $ref"
fi
