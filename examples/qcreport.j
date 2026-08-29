#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0
#
# One-pass QC report over a FASTQ file, gzipped or not.
#
#   jennifer run examples/qcreport.j reads.fastq.gz [maxReads]
#
# Inside a project that installed the deck, the import would read
# `import "@jennifer/ngs/qc.j" as qc;` - the relative paths here just let the
# example run straight from a checkout.
use io;
use os;
use convert;
import "../src/qc.j" as qc;

if (len(os.ARGS) < 2) {
    io.eprintf("usage: qcreport.j <reads.fastq[.gz]> [maxReads]\n");
    exit 2;
}

def path as string init os.ARGS[1];
def opts as qc.Options init qc.options();
if (len(os.ARGS) > 2) {
    $opts.maxReads = convert.toInt(os.ARGS[2]);
}
$opts.perBase = true;

def report as qc.Report init qc.analyze($path, $opts);
io.printf("== %s ==\n", $path);
io.printf("%s", qc.summary($report));

def means as list of float init qc.perBaseMean($report);
if (len($means) > 0) {
    io.printf("\nper-position mean quality (every 10th position)\n");
    for (def i in 0..len($means)) {
        if ($i % 10 == 0) {
            io.printf(
                "  %d|pad=4|align=right  %f|prec=1  %s\n",
                $i + 1,
                $means[$i],
                bar($means[$i]));
        }
    }
}

func bar(value as float) {
    def width as int init convert.toInt($value / 2.0);
    def out as string init "";
    for (def i in 0..$width) {
        $out = $out + "#";
    }
    return $out;
}
