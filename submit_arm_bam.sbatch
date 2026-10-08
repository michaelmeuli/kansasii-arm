#!/usr/bin/env bash
#SBATCH --job-name=arm-bam
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-130%20
#SBATCH --output=%x-%A_%a.log
# Prepare once (./run_arm_bam.sh prepare), then sbatch this array; afterwards run
#   sbatch --dependency=afterok:<arrayjob> --wrap './run_arm_bam.sh finish'  (inside env_arm)
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_arm
H=/shares/sander.imm.uzh/MM/kansasii/repos/kansasii-arm
export ARM_DIR=/shares/sander.imm.uzh/MM/kansasii/output/arm_bam
mapfile -t IDS < <(ls /shares/sander.imm.uzh/MM/kansasii/runs/mkan329/assembly/results)
"$H/06_bam_sample.sh" "${IDS[$SLURM_ARRAY_TASK_ID]}"
