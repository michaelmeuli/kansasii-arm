#!/usr/bin/env bash
# 06_bam_sample.sh - read-mapping genotyping of one sample: map reads to refs/reference.fna,
# extract the target genes (03_sample_seqs.sh bam), keep only a target-region BAM.
#
# Usage: ARM_DIR=<bam-mode dir> ./06_bam_sample.sh SAMPLE_ID [THREADS]
# Inputs : $ROOT/runs/mkan329/assembly/results/<id>/0_trimming/<id>_r{1,2}.fastq.gz
# Outputs: $ARM/work/samples/<id>/ (as 03), $ARM/bam/<id>.targets.bam (+ .bai)
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

ID=${1:?sample id}
THREADS=${2:-${SLURM_CPUS_PER_TASK:-4}}
command -v bwa >/dev/null 2>&1 || export PATH=$PATH:/home/mimeul/data/conda/envs/snippy/bin
R=$ROOT/runs/mkan329/assembly/results/$ID/0_trimming
REF=$ARM/refs/reference.fna
mkdir -p "$ARM/bam"

[[ -s "$REF.bwt" ]] || { echo "index missing: run run_arm_bam.sh first" >&2; exit 1; }
TMP=$(mktemp -d "${TMPDIR:-/tmp}/arm_bam.XXXXXX"); trap 'rm -rf "$TMP"' EXIT

bwa mem -t "$THREADS" "$REF" "$R/${ID}_r1.fastq.gz" "$R/${ID}_r2.fastq.gz" 2>/dev/null \
  | samtools sort -@ 2 -o "$TMP/$ID.bam" -
samtools index "$TMP/$ID.bam"

"$HERE/03_sample_seqs.sh" bam "$ID" "$TMP/$ID.bam"

samtools view -b -L "$ARM/work/targets.padded.bed" -o "$ARM/bam/$ID.targets.bam" "$TMP/$ID.bam"
samtools index "$ARM/bam/$ID.targets.bam"
