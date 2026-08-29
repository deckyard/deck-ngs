#!/usr/bin/env bash
#
# Guard against source headers drifting from the manifest.
#
# Every .j file in the deck opens with an SPDX pair and an interpreter floor.
# `jennifer lint` (L303) rejects a malformed pragma but says nothing about a
# missing one, so a new module can land bare and still lint clean. This checks
# that each file has all three lines, and that the two the manifest also states
# still agree with it - the licence and the engine floor are declared twice by
# design, and CI is what keeps the two copies honest.
#
# Usage: tools/check-headers.sh
set -euo pipefail

cd "$(dirname "$0")/.."
manifest="deck.yaml"

# The header block sits at the top of the file, after an optional shebang.
HEADER_LINES=5

value() {
    # value <key> - first `  key: "v"` in the manifest, unquoted.
    grep -m1 -E "^[[:space:]]+$1:" "$manifest" | cut -d: -f2- | tr -d " \"'"
}

license="$(value license)"
floor="$(value jennifer)"

if [[ -z "$license" || -z "$floor" ]]; then
    echo "error: $manifest declares no license and/or no jennifer engine" >&2
    exit 1
fi

want_license="# SPDX-License-Identifier: ${license}"
want_pragma="# pragma-jennifer-version: ${floor}"

bad=0
checked=0

while read -r file; do
    checked=$((checked + 1))
    head="$(head -n "$HEADER_LINES" "$file")"

    if ! grep -qxF -- "$want_license" <<<"$head"; then
        echo "${file}: missing '${want_license}'" >&2
        bad=$((bad + 1))
    fi
    if ! grep -qE '^# SPDX-FileCopyrightText: .' <<<"$head"; then
        echo "${file}: missing '# SPDX-FileCopyrightText: ...'" >&2
        bad=$((bad + 1))
    fi
    if ! grep -qxF -- "$want_pragma" <<<"$head"; then
        echo "${file}: missing '${want_pragma}'" >&2
        bad=$((bad + 1))
    fi
done < <(find src examples -name '*.j' | sort)

if (( bad > 0 )); then
    echo "" >&2
    echo "$bad header problem(s) across $checked files." >&2
    echo "Each .j file opens with the SPDX pair and the floor $manifest states." >&2
    exit 1
fi

echo "all $checked source files carry the SPDX pair and ${floor}"
