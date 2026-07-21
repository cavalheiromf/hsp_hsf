#!/usr/bin/env bash
#SBATCH --job-name=salmon-index
#SBATCH --partition=short
#SBATCH --array=0-3%1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --chdir=/home/mfcaval/github/hsp_hsf
#SBATCH --output=/home/mfcaval/github/hsp_hsf/logs/salmon_%A_%a.out
#SBATCH --error=/home/mfcaval/github/hsp_hsf/logs/salmon_%A_%a.err

set -euo pipefail

species=(setaria_viridis sorghum_bicolor triticum_aestivum glycine_max)
selected_species="${species[$SLURM_ARRAY_TASK_ID]}"
project_dir="${SLURM_SUBMIT_DIR:-$PWD}"
reference_dir="$project_dir/data/reference/ensembl_plants_63/$selected_species"
output_dir="$project_dir/results/rnaseq/$selected_species"
temporary_dir="$project_dir/work/salmon_transcriptomes/$selected_species"

mkdir -p "$output_dir" "$temporary_dir" "$project_dir/logs"
module load Bio/GffRead/0.12.9 Bio/Salmon/1.11.4

# gffread requires seekable, uncompressed inputs; these are temporary copies of
# the validated release-63 files and are removed after transcript extraction.
gzip -dc "$reference_dir/genome.fa.gz" > "$temporary_dir/genome.fa"
gzip -dc "$reference_dir/annotation.gff3.gz" > "$temporary_dir/annotation.gff3"

gffread "$temporary_dir/annotation.gff3" \
    -g "$temporary_dir/genome.fa" \
    -w "$output_dir/transcriptome.fa"

test -s "$output_dir/transcriptome.fa"
grep -q '^>' "$output_dir/transcriptome.fa"

salmon index \
    -t "$output_dir/transcriptome.fa" \
    -i "$output_dir/index" \
    -p "$SLURM_CPUS_PER_TASK" \
    --keepDuplicates

grep -c '^>' "$output_dir/transcriptome.fa" > "$output_dir/transcript_count.txt"
rm -f "$temporary_dir/genome.fa" "$temporary_dir/genome.fa.fai" "$temporary_dir/annotation.gff3"

echo "Completed $selected_species: $(cat "$output_dir/transcript_count.txt") transcripts"
