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

if [[ ! "$MIN_FREE_GB" =~ ^[0-9]+$ ]]; then
    echo "MIN_FREE_GB must be an integer >= 0: $MIN_FREE_GB" >&2
    exit 24
fi
if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "SLURM_CPUS_PER_TASK must be an integer > 0: $threads" >&2
    exit 25
fi
command -v flock >/dev/null || { echo "Required command not found: flock" >&2; exit 26; }
command -v mv >/dev/null || { echo "Required command not found: mv" >&2; exit 27; }

if ! row_json="$(python3 "$PROJECT_DIR/scripts/build_rnaseq_sample_manifest.py" \
  --manifest "$MANIFEST" --scope "$QUANT_SCOPE" --select-index "$task_index")"; then
    echo "Unable to select manifest row" >&2
    exit 28
fi
if ! field_text="$(python3 -c '
import json, sys
try:
    row = json.loads(sys.argv[1])
    keys = ("sample_id", "species", "r1", "r2", "salmon_index")
    values = [row[key] for key in keys]
except (json.JSONDecodeError, KeyError, TypeError) as error:
    raise SystemExit(f"Invalid manifest row: {error}")
if len(values) != 5 or not all(isinstance(value, str) and value for value in values):
    raise SystemExit("Invalid manifest row: expected five non-empty string fields")
print("\n".join(values))
' "$row_json")"; then
    echo "Unable to extract manifest fields" >&2
    exit 29
fi
mapfile -t fields <<< "$field_text"
if (( ${#fields[@]} != 5 )) || [[ -z "${fields[0]}" || -z "${fields[1]}" || -z "${fields[2]}" || -z "${fields[3]}" || -z "${fields[4]}" ]]; then
    echo "Invalid manifest row: expected exactly five non-empty fields" >&2
    exit 30
fi
sample="${fields[0]}"
species="${fields[1]}"
r1="$PROJECT_DIR/${fields[2]}"
r2="$PROJECT_DIR/${fields[3]}"
index="$PROJECT_DIR/${fields[4]}"
final="$QUANT_ROOT/$species/$sample"
partial="${final}.partial.${job_id}"
lock_file="$QUANT_ROOT/.locks/$species/$sample.lock"

echo "Started: $(date --iso-8601=seconds)"
start_epoch="$(date +%s)"
echo "Sample: $sample"
echo "Species: $species"
echo "Scope: $QUANT_SCOPE"
echo "Index: $index"

mkdir -p "$(dirname "$lock_file")"
exec {lock_fd}>"$lock_file"
flock "$lock_fd"

if [[ -e "$final" || -L "$final" ]]; then
    if [[ -d "$final" && ! -L "$final" ]] && python3 "$PROJECT_DIR/scripts/validate_salmon_quant.py" --single "$final" >/dev/null; then
        echo "$sample already complete; skipping"
        exit 0
    fi
    echo "Existing final output is invalid and will not be overwritten: $final" >&2
    exit 20
fi
if compgen -G "${final}.partial.*" > /dev/null; then
    echo "Partial output already exists and requires review: ${final}.partial.*" >&2
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
mv -T -- "$partial" "$final"
echo "Completed: $(date --iso-8601=seconds)"
