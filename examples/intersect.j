#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0
#
# bedtools-style interval arithmetic: how much of each BED feature is covered by
# the exons of a GFF/GTF annotation.
#
#   jennifer run examples/intersect.j peaks.bed genes.gff3
use io;
use os;
import "../src/bed.j" as bed;
import "../src/gff.j" as gff;
import "../src/intervals.j" as intervals;

if (len(os.ARGS) < 3) {
    io.eprintf("usage: intersect.j <peaks.bed> <annotation.gff3|gtf>\n");
    exit 2;
}

def peaks as list of intervals.Interval init bed.readAll(os.ARGS[1]);
def exons as list of intervals.Interval init gff.readIntervals(os.ARGS[2], "exon");

io.printf("peaks %d, exons %d\n", len($peaks), len($exons));
io.printf("exonic bases (merged) %d\n", intervals.totalLength($exons));

def sorted as list of intervals.Interval init intervals.sort($peaks);
def hits as list of int init intervals.countOverlaps($peaks, $exons);
def covered as list of int init intervals.coverage($peaks, $exons);

io.printf("\n%s\t%s\t%s\t%s\n", "peak", "exons", "coveredBp", "fraction");
for (def i in 0..len($sorted)) {
    def iv as intervals.Interval init $sorted[$i];
    def width as int init intervals.length($iv);
    def fraction as float init 0.0;
    if ($width > 0) {
        $fraction = $covered[$i] / $width;
    }
    io.printf("%s\t%d\t%d\t%f|prec=3\n", intervals.toBed3($iv), $hits[$i], $covered[$i], $fraction);
}

def shared as list of intervals.Interval init intervals.intersect($peaks, $exons);
io.printf("\noverlap regions %d, total %d bp\n", len($shared), intervals.totalLength($shared));

def outside as list of intervals.Interval init intervals.subtract($peaks, $exons);
io.printf(
    "non-exonic peak pieces %d, total %d bp\n",
    len($outside),
    intervals.totalLength($outside));
