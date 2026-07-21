#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 SPECIES_ID FAMILY [CPUS]" >&2
    exit 2
fi

species="$1"
family="$2"
cpus="${3:-4}"
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
families_config="$project_dir/config/pfam_families.tsv"
proteome="$project_dir/work/representative_proteomes/${species}.fa"
models_dir="$project_dir/data/pfam/38.2/models"
outdir="$project_dir/work/hmmer/$species"

pfam_id="$(awk -F '\t' -v family="$family" 'NR > 1 && $1 == family {print $2}' "$families_config")"
if [[ -z "$pfam_id" ]]; then
    echo "Unknown family: $family" >&2
    exit 3
fi

model="$models_dir/${family}__${pfam_id}.hmm"
if [[ ! -s "$proteome" || ! -s "$model" ]]; then
    echo "Missing proteome or HMM: $proteome | $model" >&2
    exit 4
fi

if ! command -v hmmsearch >/dev/null 2>&1; then
    module load Bio/HMMER3/3.4
fi

mkdir -p "$outdir"
hmmsearch \
    --cut_ga \
    --cpu "$cpus" \
    --domtblout "$outdir/${family}__${pfam_id}.domtblout" \
    "$model" \
    "$proteome" \
    > "$outdir/${family}__${pfam_id}.report.txt"

test -s "$outdir/${family}__${pfam_id}.domtblout"
