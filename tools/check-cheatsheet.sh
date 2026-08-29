#!/usr/bin/env bash
#
# Guard against docs/cheatsheet.md drifting from the code.
#
# The cheatsheet is hand-committed rather than generated at build time, so
# this checks the one property that matters: every exported name in every module
# appears in it. Adding an export without documenting it fails CI.
#
# Usage: tools/check-cheatsheet.sh
set -euo pipefail

cd "$(dirname "$0")/.."
ref="docs/cheatsheet.md"

if [[ ! -f "$ref" ]]; then
    echo "error: $ref not found" >&2
    exit 1
fi

missing=0
checked=0

for file in src/*.j; do
    case "$file" in
        *_test.j) continue ;;
    esac
    mod="$(basename "${file}" .j)"

    while read -r kind name; do
        checked=$((checked + 1))
        if [[ "$kind" == "func" ]]; then
            token="\`${mod}.${name}("
        else
            token="\`${mod}.${name}\`"
        fi
        if ! grep -qF -- "$token" "$ref"; then
            echo "undocumented export: ${mod}.${name} (${kind})" >&2
            missing=$((missing + 1))
        fi
    done < <(
        grep -oE '^export (func|def const|def struct|def enum) [A-Za-z][A-Za-z0-9_]*' "$file" \
            | sed -E 's/^export (func|def const|def struct|def enum) /\1 /; s/^def //'
    )
done

if (( missing > 0 )); then
    echo "" >&2
    echo "$missing of $checked exported names are missing from $ref." >&2
    echo "Add them there, then re-run tools/check-cheatsheet.sh." >&2
    exit 1
fi

echo "all $checked exported names are documented in $ref"
