#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0
#
# Stream a SAM file, count alignment categories and write primary alignments
# above a mapping-quality cutoff out as BED.
#
#   jennifer run examples/samfilter.j alignments.sam 30 out.bed
#   samtools view -h aln.bam | jennifer run examples/samfilter.j - 30 out.bed
use io;
use os;
use convert;
import "../src/sam.j" as sam;
import "../src/bed.j" as bed;
import "../src/source.j" as source;
import "../src/intervals.j" as intervals;

if (len(os.ARGS) < 4) {
    io.eprintf("usage: samfilter.j <alignments.sam|-> <minMapq> <out.bed>\n");
    exit 2;
}

def minMapq as int init convert.toInt(os.ARGS[2]);
def reader as sam.Reader init openInput(os.ARGS[1]);
def sink as source.Sink init source.create(os.ARGS[3]);

def total as int init 0;
def unmapped as int init 0;
def secondary as int init 0;
def duplicate as int init 0;
def lowMapq as int init 0;
def written as int init 0;

def going as bool init true;
while ($going) {
    def step as sam.Next init sam.next($reader);
    match ($step) {
        when Alignment(hit) { $going = tally($hit.value); }
        when Done { $going = false; }
    }
}
sam.close($reader);
def code as int init source.closeSink($sink);

io.printf("alignments   %d\n", $total);
io.printf("unmapped     %d\n", $unmapped);
io.printf("secondary    %d\n", $secondary);
io.printf("duplicate    %d\n", $duplicate);
io.printf("mapq < %d    %d\n", $minMapq, $lowMapq);
io.printf("written      %d\n", $written);

func openInput(path as string) {
    if ($path == "-") {
        return sam.openStdin();
    }
    return sam.open($path);
}

func tally(rec as sam.Record) {
    $total = $total + 1;
    if (sam.isUnmapped($rec)) {
        $unmapped = $unmapped + 1;
        return true;
    }
    if (sam.isSecondary($rec) or sam.isSupplementary($rec)) {
        $secondary = $secondary + 1;
        return true;
    }
    if (sam.isDuplicate($rec)) {
        $duplicate = $duplicate + 1;
        return true;
    }
    if ($rec.mapq < $minMapq) {
        $lowMapq = $lowMapq + 1;
        return true;
    }
    def iv as intervals.Interval init sam.toInterval($rec);
    bed.write($sink, $iv);
    $written = $written + 1;
    return true;
}
