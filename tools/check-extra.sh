#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
#
# The repo-specific checks. Both workflows call this if it is executable.
#
# ngs is a modules deck, so tools/inject-version.sh (topics decks only) does
# not stamp its VERSION: src/ngs.j keeps it by hand. Releases 0.2.0 and 0.3.0
# shipped a stale constant ("0.1.0", "0.3.0-dev") because only the manifest
# was bumped. Fail while the two differ.
#
# Usage: tools/check-extra.sh
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck source=tools/deck-lib.sh
source tools/deck-lib.sh
deck_load

constant="$(sed -n 's/^export def const VERSION as string init "\(.*\)";$/\1/p' src/ngs.j)"
if [[ -z "$constant" ]]; then
    deck_die "src/ngs.j has no VERSION constant"
fi
if [[ "$constant" != "$DECK_VERSION" ]]; then
    deck_die "src/ngs.j VERSION is \"${constant}\" but ${DECK_MANIFEST} declares \"${DECK_VERSION}\"; bump them together"
fi
echo "src/ngs.j VERSION matches ${DECK_MANIFEST} (${DECK_VERSION})"
