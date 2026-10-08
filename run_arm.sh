#!/usr/bin/env bash
# run_arm.sh - run the full AMR-target pipeline on the Mkan329 assemblies.
#
# Usage: ./run_arm.sh [THREADS]
#
# Inputs : $KANSASII_ROOT/runs/mkan329/assembly/results/<id>/1_unicycler/<id>.fasta
# Outputs: $ARM_DIR (default $KANSASII_ROOT/output/arm)/{refs,work}
# Env    : MIN_DEPTH / MIN_IDENT as in 03_sample_seqs.sh; ASM_GLOB overrides the input glob.
# Samples already present in work/samples/ are skipped (delete the dir to redo one);
# 04/05 always re-run on everything present.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

THREADS=${1:-${SLURM_CPUS_PER_TASK:-4}}
ASM_GLOB=${ASM_GLOB:-$ROOT/runs/mkan329/assembly/results/*/1_unicycler/*.fasta}

[[ -s "$ARM/refs/reference.gff" ]] || "$HERE/01_fetch_reference.sh"
python3 "$HERE/02_make_targets.py"

n=0
for fa in $ASM_GLOB; do
  [[ -e "$fa" ]] || { echo "no assemblies match: $ASM_GLOB" >&2; exit 1; }
  id=$(basename "$fa" .fasta)
  if [[ -d "$ARM/work/samples/$id/genes" ]]; then echo "[skip] $id"; continue; fi
  "$HERE/03_sample_seqs.sh" asm "$id" "$fa"
  n=$((n+1))
done
echo "processed $n new samples"

"$HERE/04_align.sh" "$THREADS"
python3 "$HERE/05_annotate.py"
