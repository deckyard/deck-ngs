#!/usr/bin/env -S jennifer run
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0
#
# Stream a VCF, keep PASS SNVs above a quality cutoff, and report per-sample
# genotype counts.
#
#   jennifer run examples/vcffilter.j calls.vcf.gz 30
use io;
use os;
use convert;
use maps;
import "../src/vcf.j" as vcf;

if (len(os.ARGS) < 3) {
    io.eprintf("usage: vcffilter.j <calls.vcf[.gz]> <minQual>\n");
    exit 2;
}

def minQual as float init convert.toFloat(os.ARGS[2]);
def reader as vcf.Reader init vcf.open(os.ARGS[1]);

def total as int init 0;
def kept as int init 0;
def snvs as int init 0;
def indels as int init 0;
def genotypes as map of string to int init {};

def going as bool init true;
while ($going) {
    def step as vcf.Next init vcf.next($reader);
    match ($step) {
        when Record(hit) { $going = consider($hit.value); }
        when Done { $going = false; }
    }
}
vcf.closeStrict($reader);

io.printf("samples      %s\n", joinNames($reader.samples));
io.printf("variants     %d\n", $total);
io.printf("SNVs         %d\n", $snvs);
io.printf("indels       %d\n", $indels);
io.printf("kept         %d\n", $kept);
for (def call in $genotypes) {
    io.printf("  %s %d\n", $call, $genotypes[$call]);
}

func consider(v as vcf.Variant) {
    $total = $total + 1;
    if (vcf.isSnv($v)) {
        $snvs = $snvs + 1;
    }
    if (vcf.isIndel($v)) {
        $indels = $indels + 1;
    }
    if (not vcf.passed($v) or $v.qual < $minQual or not vcf.isSnv($v)) {
        return true;
    }
    $kept = $kept + 1;
    for (def i in 0..len($reader.samples)) {
        def call as string init vcf.genotype($v, $i);
        if (maps.has($genotypes, $call)) {
            $genotypes[$call] = $genotypes[$call] + 1;
        } else {
            $genotypes[$call] = 1;
        }
    }
    return true;
}

func joinNames(names as list of string) {
    if (len($names) == 0) {
        return "(none)";
    }
    def out as string init $names[0];
    for (def i in 1..len($names)) {
        $out = $out + ", " + $names[$i];
    }
    return $out;
}
