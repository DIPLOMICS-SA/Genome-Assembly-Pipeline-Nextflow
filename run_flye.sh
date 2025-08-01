#!/bin/bash
set -euo pipefail

# Args passed from master.sh
READS=${1:?Provide FASTQ}
GENOME_SIZE=${2:?Provide genome size}
COVERAGE=${3:?Provide asm coverage}
THREADS=${4:-16}
READ_TYPE=${5:-nano-raw}

module purge
module load chpc/BIOMODULES
module load flye/2.9.5

# 🔧 Fix the relative path issue
READS=$(realpath "$READS")

OUTDIR="results/Flye_results"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

if [[ -f params.json && -d 00-assembly ]]; then
    echo "🔁 Resuming Flye..."
    flye --${READ_TYPE} "$READS" -o . \
         --resume \
         --genome-size "$GENOME_SIZE" \
         --asm-coverage "$COVERAGE" \
         -t "$THREADS"
else
    echo "🚀 Running Flye from scratch..."
    flye --${READ_TYPE} "$READS" -o . \
         --genome-size "$GENOME_SIZE" \
         --asm-coverage "$COVERAGE" \
         -t "$THREADS"
fi

