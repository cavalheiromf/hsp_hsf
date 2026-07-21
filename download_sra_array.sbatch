#!/usr/bin/env bash
#SBATCH --job-name=sra_download
#SBATCH --partition=long
#SBATCH --array=1-58%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/home/mfcaval/github/hsp_hsf
#SBATCH --output=/home/mfcaval/github/hsp_hsf/logs/sra_%A_%a.out
#SBATCH --error=/home/mfcaval/github/hsp_hsf/logs/sra_%A_%a.err

set -euo pipefail

# Each array task handles one biological run from the curated SRR list.
PROJECT_DIR="/home/mfcaval/github/hsp_hsf"
SRR_LIST="$PROJECT_DIR/reference/srr_list.txt"
SRA_DIR="$PROJECT_DIR/data/sra"
FASTQ_DIR="$PROJECT_DIR/data/fastq"
TMP_ROOT="$PROJECT_DIR/data/tmp"
SRA_BIN="/home/mfcaval/.local/share/mamba/envs/sraEnv/bin"

PREFETCH="$SRA_BIN/prefetch"
VDB_VALIDATE="$SRA_BIN/vdb-validate"
FASTERQ_DUMP="$SRA_BIN/fasterq-dump"
PIGZ="/usr/bin/pigz"

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMP_ROOT"

for executable in "$PREFETCH" "$VDB_VALIDATE" "$FASTERQ_DUMP" "$PIGZ"; do
    if [[ ! -x "$executable" ]]; then
        echo "Required executable not found: $executable" >&2
        exit 10
    fi
done

# SLURM_ARRAY_TASK_ID is the one-based line number in srr_list.txt.
srr="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")"

if [[ ! "$srr" =~ ^SRR[0-9]+$ ]]; then
    echo "Invalid SRR on line ${SLURM_ARRAY_TASK_ID}: '$srr'" >&2
    exit 11
fi

task_tmp="${SLURM_TMPDIR:-${TMP_ROOT}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$task_tmp"

sra_file="$SRA_DIR/$srr/$srr.sra"
r1="$FASTQ_DIR/${srr}_1.fastq"
r2="$FASTQ_DIR/${srr}_2.fastq"
r1_gz="${r1}.gz"
r2_gz="${r2}.gz"

echo "Started:     $(date --iso-8601=seconds)"
echo "Job:         $SLURM_JOB_ID"
echo "Array task:  $SLURM_ARRAY_TASK_ID"
echo "Node:        $(hostname)"
echo "SRR:         $srr"
echo "Temporary:   $task_tmp"

# A completed pair is idempotent: resubmitting the array will skip it.
if [[ -s "$r1_gz" && -s "$r2_gz" ]]; then
    echo "$srr already has both compressed FASTQ files; skipping."
    exit 0
fi

# Do not silently overwrite a partially converted sample.
if [[ -e "$r1" || -e "$r2" || -e "$r1_gz" || -e "$r2_gz" ]]; then
    echo "Partial FASTQ output exists for $srr; manual review required." >&2
    exit 12
fi

# Prefetch can resume an interrupted accession download in the same directory.
"$PREFETCH" \
    --max-size 200g \
    --output-directory "$SRA_DIR" \
    "$srr"

"$VDB_VALIDATE" "$sra_file"

# Split the paired-end RNA-seq run into R1 and R2 FASTQ files.
"$FASTERQ_DUMP" \
    "$sra_file" \
    --split-files \
    --threads "$SLURM_CPUS_PER_TASK" \
    --temp "$task_tmp" \
    --outdir "$FASTQ_DIR"

if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    echo "FASTQ conversion did not produce a non-empty R1/R2 pair for $srr." >&2
    exit 13
fi

"$PIGZ" -p "$SLURM_CPUS_PER_TASK" "$r1" "$r2"

if [[ ! -s "$r1_gz" || ! -s "$r2_gz" ]]; then
    echo "Compression did not produce a non-empty R1/R2 pair for $srr." >&2
    exit 14
fi

echo "Completed:   $(date --iso-8601=seconds)"
ls -lh "$r1_gz" "$r2_gz"
