#!/usr/bin/env bash
#SBATCH --job-name=hsp_hmmer
#SBATCH --partition=short
#SBATCH --array=0-27%4
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/mfcaval/github/hsp_hsf
#SBATCH --output=/home/mfcaval/github/hsp_hsf/logs/hmmer/hmmer_%A_%a.out
#SBATCH --error=/home/mfcaval/github/hsp_hsf/logs/hmmer/hmmer_%A_%a.err

set -euo pipefail

species=(setaria_viridis sorghum_bicolor triticum_aestivum glycine_max)
families=(Hsp20 Hsp40 Hsp60 Hsp70 Hsp90 Hsp100 HSF)

species_index=$((SLURM_ARRAY_TASK_ID / ${#families[@]}))
family_index=$((SLURM_ARRAY_TASK_ID % ${#families[@]}))

selected_species="${species[$species_index]}"
selected_family="${families[$family_index]}"

echo "Started: $(date --iso-8601=seconds)"
echo "Species: $selected_species"
echo "Family:  $selected_family"

scripts/run_hmmer.sh \
    "$selected_species" \
    "$selected_family" \
    "$SLURM_CPUS_PER_TASK"

echo "Completed: $(date --iso-8601=seconds)"
