# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0
/**
 * White-box tests for the deck's entry module.
 * Run with: jennifer test ngs_test.j
 */
use testing;
use strings;

func testIdentityMatchesTheManifest() {
    testing.assertEqual(NAME, "@deckyard/ngs");
}

func testCapabilitiesAreConsistent() {
    def caps as Capabilities init capabilities();
    # Streaming gzip is external tooling plus a gzip on PATH; it can never be
    # reported as available when either half is missing.
    if ($caps.streamingGzip) {
        testing.assertEqual($caps.externalTools, true);
        testing.assertEqual($caps.gzipOnPath, true);
    }
}

func testDescribeNamesTheDeck() {
    def text as string init describe();
    testing.assertEqual(strings.startsWith($text, NAME + " ("), true);
}

func testRequireFullHostMatchesCapabilities() {
    if (capabilities().externalTools) {
        requireFullHost();
        return;
    }
    testing.assertThrows("requireFullHost", "ngs.ngs");
}
