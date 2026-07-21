#!/usr/bin/env bash
#SBATCH --job-name=salmon-quant
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --chdir=/home/mfcaval/github/hsp_hsf
#SBATCH --output=/home/mfcaval/github/hsp_hsf/logs/salmon_quant_%A_%a.out
#SBATCH --error=/home/mfcaval/github/hsp_hsf/logs/salmon_quant_%A_%a.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/mfcaval/github/hsp_hsf}"
MANIFEST="${MANIFEST:-$PROJECT_DIR/config/rnaseq_samples.tsv}"
QUANT_ROOT="${QUANT_ROOT:-$PROJECT_DIR/results/rnaseq/quant}"
QUANT_SCOPE="${QUANT_SCOPE:-canary}"
MIN_FREE_GB="${MIN_FREE_GB:-5}"
task_index="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
threads="${SLURM_CPUS_PER_TASK:-16}"
job_id="${SLURM_JOB_ID:-manual}"

row_json="$(python3 "$PROJECT_DIR/scripts/build_rnaseq_sample_manifest.py" \
  --manifest "$MANIFEST" --scope "$QUANT_SCOPE" --select-index "$task_index")"
mapfile -t fields < <(python3 -c '
import json, sys
row = json.loads(sys.argv[1])
for key in ("sample_id", "species", "r1", "r2", "salmon_index"):
    print(row[key])
' "$row_json")
sample="${fields[0]}"
species="${fields[1]}"
r1="$PROJECT_DIR/${fields[2]}"
r2="$PROJECT_DIR/${fields[3]}"
index="$PROJECT_DIR/${fields[4]}"
final="$QUANT_ROOT/$species/$sample"
partial="${final}.partial.${job_id}"

echo "Started: $(date --iso-8601=seconds)"
start_epoch="$(date +%s)"
echo "Sample: $sample"
echo "Species: $species"
echo "Scope: $QUANT_SCOPE"
echo "Index: $index"

if [[ -d "$final" ]]; then
    if python3 "$PROJECT_DIR/scripts/validate_salmon_quant.py" --single "$final" >/dev/null; then
        echo "$sample already complete; skipping"
        exit 0
    fi
    echo "Existing final output is invalid and will not be overwritten: $final" >&2
    exit 20
fi
if [[ -e "$partial" ]]; then
    echo "Partial output already exists and requires review: $partial" >&2
    exit 21
fi
for path in "$r1" "$r2" "$index/versionInfo.json"; do
    [[ -s "$path" ]] || { echo "Missing or empty input: $path" >&2; exit 22; }
done
gzip -t "$r1"
gzip -t "$r2"
mkdir -p "$(dirname "$final")"
free_kb="$(df -Pk "$(dirname "$final")" | awk 'NR==2 {print $4}')"
required_kb="$((MIN_FREE_GB * 1024 * 1024))"
if (( free_kb < required_kb )); then
    echo "Insufficient free space: ${free_kb} KiB available, ${required_kb} KiB required" >&2
    exit 23
fi

if [[ -z "${SALMON_BIN:-}" ]]; then
    module load Bio/Salmon/1.11.4
    SALMON_BIN="$(command -v salmon)"
fi
command=(
    "$SALMON_BIN" quant
    -i "$index"
    -l A
    -1 "$r1"
    -2 "$r2"
    --validateMappings
    --seqBias
    --gcBias
    -p "$threads"
    -o "$partial"
)
printf 'Command:'
printf ' %q' "${command[@]}"
printf '\n'
"${command[@]}"
end_epoch="$(date +%s)"
printf '{"elapsed_seconds": %d}\n' "$((end_epoch - start_epoch))" > "$partial/run_metrics.json"
python3 "$PROJECT_DIR/scripts/validate_salmon_quant.py" --single "$partial"
mv "$partial" "$final"
echo "Completed: $(date --iso-8601=seconds)"
