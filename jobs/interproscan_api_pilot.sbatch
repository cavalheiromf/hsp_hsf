#!/usr/bin/env bash
#SBATCH --job-name=ipr-api
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=24:00:00
#SBATCH --output=logs/interproscan_api_%A_%a.out
#SBATCH --error=logs/interproscan_api_%A_%a.err
#SBATCH --array=0-3%1

set -euo pipefail

project_dir="${SLURM_SUBMIT_DIR:-$PWD}"
species=(setaria_viridis sorghum_bicolor triticum_aestivum glycine_max)
selected_species="${species[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "${INTERPRO_EMAIL:-}" ]]; then
    echo "Set INTERPRO_EMAIL before submitting this job." >&2
    exit 2
fi

mkdir -p "$project_dir/logs" "$project_dir/results/interproscan"

python3 "$project_dir/scripts/run_interproscan_api.py" \
    --input "$project_dir/work/candidates/${selected_species}.fa" \
    --species "$selected_species" \
    --outdir "$project_dir/results/interproscan" \
    --email "$INTERPRO_EMAIL" \
    --chunk-size 100 \
    --poll-seconds 30
