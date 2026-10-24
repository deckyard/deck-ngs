# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0
/**
 * White-box tests for the deck's entry module.
 * Run with: jennifer test ngs_test.j
 */
use testing;
use strings;

func testIdentityMatchesTheManifest() {
    testing.assertEqual(NAME, "@deckyard/ngs");
    testing.assertEqual(VERSION, "0.1.0");
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

func testDescribeNamesTheDeckAndVersion() {
    def text as string init describe();
    testing.assertEqual(strings.contains($text, NAME), true);
    testing.assertEqual(strings.contains($text, VERSION), true);
}

func testRequireFullHostMatchesCapabilities() {
    if (capabilities().externalTools) {
        requireFullHost();
        return;
    }
    testing.assertThrows("requireFullHost", "ngs.ngs");
}
