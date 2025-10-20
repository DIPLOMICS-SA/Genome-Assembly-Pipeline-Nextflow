#!/bin/bash
set -euo pipefail

# ==============================
#   Flye Assembly Script
#   (Containerized version)
# ==============================

READS=${1:?Provide FASTQ}
GENOME_SIZE=${2:?Provide genome size}
COVERAGE=${3:?Provide asm coverage}
THREADS=${4:-16}
READ_TYPE=${5:-nano-raw}
SPECIES_NAME=${6:?Provide species name}

CONTAINER="/home/apps/chpc/bio/1ksa_pipeline/1ksa_pipeline.sif"

module purge
module load chpc/BIOMODULES
module load apptainer/1.2.3_SUID

READS=$(realpath "$READS")

OUTDIR=$(realpath "results/Flye_results")
mkdir -p "$OUTDIR"
cd "$OUTDIR"

LOG="${SPECIES_NAME}_flye.log"

if [[ -f params.json && -d 00-assembly ]]; then
    echo "Resuming Flye assembly..."
    apptainer exec --bind /mnt/lustre,/home "$CONTAINER" \
        flye --${READ_TYPE} "$READS" -o . \
             --resume \
             --genome-size "$GENOME_SIZE" \
             --asm-coverage "$COVERAGE" \
             -t "$THREADS" > "$LOG" 2>&1
else
    echo "Running Flye assembly from scratch..."
    apptainer exec --bind /mnt/lustre,/home "$CONTAINER" \
        flye --${READ_TYPE} "$READS" -o . \
             --genome-size "$GENOME_SIZE" \
             --asm-coverage "$COVERAGE" \
             -t "$THREADS" > "$LOG" 2>&1
fi

if [[ -f assembly.fasta ]]; then
    mv assembly.fasta "${SPECIES_NAME}_assembly.fasta"
fi
