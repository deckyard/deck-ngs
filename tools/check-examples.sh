#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
#
# ngs overrides the shared example runner.
#
# Every example takes input files - a FASTQ, a BED, a GFF3, a SAM, a VCF - so
# there is nothing to run without them. This builds a tiny fixture set in a
# temporary directory, runs each example against it, and cleans up. The point
# is that the examples still work end to end, not that they produce any
# particular numbers.
#
# examples/pipeline.j drives external aligners, which are not on a runner, so
# it is parse-checked rather than run.
#
# The fixtures live *inside the repository*, and every path handed to an
# example is *relative* to it. Both matter: in CI `jennifer` is a shim around a
# container that mounts only the workspace, at a different absolute path than
# the runner sees. A fixture under /tmp is invisible in there, and an absolute
# host path resolves to nothing. Relative paths under the checkout are the only
# form that works both natively and through the container.
#
# Usage: tools/check-examples.sh
set -euo pipefail
cd "$(dirname "$0")/.."

JENNIFER="${JENNIFER:-jennifer}"
work=".ci-fixtures"
rm -rf "$work"
mkdir -p "$work"
trap 'rm -rf "$work"' EXIT

# --- fixtures --------------------------------------------------------------
{
    printf "@r1 first\nACGTACGTNN\n+\nIIIIIIII##\n"
    printf "@r2 second\nGGCCGGCCAT\n+\nIIIIIIIIII\n"
    printf "@r3 third\nTTTTAAAACG\n+\n##########\n"
} > "$work/reads.fastq"
gzip -kf "$work/reads.fastq"

printf "chr1\t100\t200\tpeak1\t55\t+\nchr1\t400\t500\tpeak2\t20\t-\n" > "$work/peaks.bed"
printf "##gff-version 3\nchr1\t.\tgene\t150\t450\t.\t+\t.\tID=g1\nchr1\t.\texon\t150\t250\t.\t+\t.\tID=e1\n" > "$work/genes.gff3"
{
    printf "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
    printf "r1\t0\tchr1\t101\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
    printf "r2\t16\tchr1\t401\t20\t10M\t*\t0\t0\tGGCCGGCCAT\tIIIIIIIIII\n"
} > "$work/aln.sam"
{
    printf "##fileformat=VCFv4.2\n"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    printf "chr1\t150\t.\tA\tG\t60\tPASS\tDP=30\tGT\t0/1\n"
    printf "chr1\t420\t.\tACGT\tA\t10\tPASS\tDP=4\tGT\t1/1\n"
} > "$work/calls.vcf"

# --- one invocation per example --------------------------------------------
run() {
    echo "--- $1"
    shift
    "$JENNIFER" run "$@" > /dev/null
}

failed=0
run examples/qcreport.j  examples/qcreport.j  "$work/reads.fastq.gz" 1000        || failed=1
run examples/trimreads.j examples/trimreads.j "$work/reads.fastq.gz" "$work/out.fastq.gz" || failed=1
run examples/intersect.j examples/intersect.j "$work/peaks.bed" "$work/genes.gff3"        || failed=1
run examples/samfilter.j examples/samfilter.j "$work/aln.sam" 30 "$work/hits.bed"         || failed=1
run examples/vcffilter.j examples/vcffilter.j "$work/calls.vcf" 30                        || failed=1

# Drives bwa/samtools, which a runner does not have: parse it instead.
echo "--- examples/pipeline.j (parse only - needs external aligners)"
"$JENNIFER" ast examples/pipeline.j > /dev/null || failed=1

if (( failed )); then
    echo "" >&2
    echo "at least one example failed." >&2
    exit 1
fi
echo "all 6 examples ran"
