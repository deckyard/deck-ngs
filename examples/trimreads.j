#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0
#
# Adapter- and quality-trim a FASTQ file, streaming in and out.
#
#   jennifer run examples/trimreads.j in.fastq.gz out.fastq.gz
#
# The output is gzip-compressed automatically because the name ends in `.gz`.
use io;
use os;
import "../src/trim.j" as trim;

if (len(os.ARGS) < 3) {
    io.eprintf("usage: trimreads.j <in.fastq[.gz]> <out.fastq[.gz]>\n");
    exit 2;
}

def opts as trim.Options init trim.options();
$opts.qualityThreshold = 20;
$opts.minLength = 36;
$opts.maxAmbiguous = 5;

def st as trim.Stats init trim.run(os.ARGS[1], os.ARGS[2], $opts);

io.printf("reads in         %d\n", $st.readsIn);
io.printf("reads out        %d\n", $st.readsOut);
io.printf("bases in         %d\n", $st.basesIn);
io.printf("bases out        %d\n", $st.basesOut);
io.printf("adapter trimmed  %d\n", $st.adapterTrimmed);
io.printf("quality trimmed  %d\n", $st.qualityTrimmed);
io.printf(
    "dropped: short %d, ambiguous %d, low quality %d\n",
    $st.tooShort,
    $st.tooAmbiguous,
    $st.lowQuality);
