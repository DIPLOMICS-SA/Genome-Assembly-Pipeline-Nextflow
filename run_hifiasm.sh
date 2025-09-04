#!/bin/bash
set -euo pipefail

SPECIES_NAME=${2:?Provide species name}
THREADS=${3:-16}

module purge
module load chpc/BIOMODULES
module load hifiasm/0.24.0-r703

# Path to trimmed FASTQ
TRIMMED=$(realpath "results/trimmed_fastq/${SPECIES_NAME}.trimmed.fastq")

if [[ ! -f "$TRIMMED" ]]; then
    echo "  ERROR: Trimmed FASTQ file not found at $TRIMMED"
    exit 1
fi

OUTDIR=$(realpath "results/Hifiasm_results")
mkdir -p "$OUTDIR"
cd "$OUTDIR"

LOG="${SPECIES_NAME}_hifiasm.log"

if [[ -f ${SPECIES_NAME}.bp.p_ctg.gfa ]]; then
    echo "✅ Hifiasm assembly already exists, skipping."
    exit 0
fi

echo "🚀 Running Hifiasm..."
hifiasm -t "$THREADS" --ont "$TRIMMED" -o "$SPECIES_NAME" > "$LOG" 2>&1

# Convert GFA to FASTA with species name prefix
for file in ${SPECIES_NAME}.bp.p_ctg.gfa ${SPECIES_NAME}.bp.hap1.p_ctg.gfa ${SPECIES_NAME}.bp.hap2.p_ctg.gfa; do
    if [[ -f "$file" ]]; then
        suffix=$(echo "$file" | cut -d'.' -f3 | sed 's/^p_//')
        output="${SPECIES_NAME}_${suffix}.fasta"
        awk '/^S/{print ">"$2; print $3}' "$file" > "$output"
    fi
done
