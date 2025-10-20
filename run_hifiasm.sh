#!/bin/bash
set -euo pipefail

# ==============================
#   Hifiasm Assembly Script
#   (Containerized version)
# ==============================

SPECIES_NAME=${1:?Provide species name}
THREADS=${2:-16}

# Path to container
CONTAINER="/home/apps/chpc/bio/1ksa_pipeline/1ksa_pipeline.sif"

# Load Apptainer
module purge
module load chpc/BIOMODULES
module load apptainer/1.2.3_SUID

# Path to trimmed FASTQ
TRIMMED=$(realpath "results/trimmed_fastq/${SPECIES_NAME}.trimmed.fastq")

if [[ ! -f "$TRIMMED" ]]; then
    echo "ERROR: Trimmed FASTQ file not found at $TRIMMED"
    exit 1
fi

# Create output directory
OUTDIR=$(realpath "results/Hifiasm_results")
mkdir -p "$OUTDIR"
cd "$OUTDIR"

LOG="${SPECIES_NAME}_hifiasm.log"

# If assembly already exists, skip
if [[ -f ${SPECIES_NAME}.bp.p_ctg.gfa ]]; then
    echo "Hifiasm assembly already exists, skipping."
    exit 0
fi

# Run hifiasm inside container
echo "Running Hifiasm..."
apptainer exec --bind /mnt/lustre,/home "$CONTAINER" \
    hifiasm -t "$THREADS" --ont "$TRIMMED" -o "$SPECIES_NAME" > "$LOG" 2>&1

# Convert GFA outputs to FASTA with species name prefixes
for file in ${SPECIES_NAME}.bp.p_ctg.gfa ${SPECIES_NAME}.bp.hap1.p_ctg.gfa ${SPECIES_NAME}.bp.hap2.p_ctg.gfa; do
    if [[ -f "$file" ]]; then
        suffix=$(echo "$file" | cut -d'.' -f3 | sed 's/^p_//')
        output="${SPECIES_NAME}_${suffix}.fasta"
        echo "📄 Converting $file → $output"
        awk '/^S/{print ">"$2; print $3}' "$file" > "$output"
    fi
done

echo "Hifiasm assembly complete for ${SPECIES_NAME}"