#!/bin/bash

module add chpc/BIOMODULES
module add blobtools

# Usage: ./blobtools_snail.sh species_name

SPECIES=$1

if [ -z "$SPECIES" ]; then
    echo "Usage: $0 <species_name>"
    exit 1
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="./results"
RESULTS_DIR="${BASE}/${SPECIES}"
BUSCO_DIR="${BASE}/${SPECIES}_other_results_outputs/Busco_results/Busco_output"

# Locate files
FASTA=$(find "$RESULTS_DIR" -name "*_racon_polished.fasta" | head -n 1)
BUSCO=$(find "$BUSCO_DIR" -name "full_table.tsv" | head -n 1)

# Check files exist
[[ -f "$FASTA" ]] || { echo "ERROR: Racon polished FASTA not found."; exit 1; }
[[ -f "$BUSCO" ]] || { echo "ERROR: BUSCO full_table.tsv not found."; exit 1; }

# Create BlobDir
OUTDIR="${SPECIES}_blobtools"
mkdir -p "$OUTDIR"

cp "$FASTA" "$OUTDIR/"
cp "$BUSCO" "$OUTDIR/"

singularity run "$SIF" blobtools add \
    --fasta "$OUTDIR/$(basename "$FASTA")" \
    --busco "$OUTDIR/$(basename "$BUSCO")" \
    --create \
    "$OUTDIR"

singularity run "$SIF" blobtools view \
    --plot \
    --view snail \
    "$OUTDIR"

echo "Done! Snail plot created in ${OUTDIR}"
