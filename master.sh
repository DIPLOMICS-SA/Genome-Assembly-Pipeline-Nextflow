#!/bin/bash
set -euo pipefail

CONFIG="params.config"
if [[ ! -f $CONFIG ]]; then
    echo "Missing $CONFIG file."
    exit 1
fi
source "$CONFIG"

# Ensure assembler is set in params.config
if [ -z "${assembler:-}" ]; then
    echo "ERROR: 'assembler' is not set in params.config"
    exit 1
fi

FASTQ=${1:?Please provide FASTQ input file}
ASSEMBLER="${2:-$assembler}"
LINEAGE="${LINEAGE:?LINEAGE not set in params.config}"

# Use species_name from params.config
TRIMMED="results/trimmed_fastq/${species_name}.trimmed.fastq"
QC1="results/nanoplot_before_trim/NanoPlot_CHECK_1/NanoStats.txt"
QC2="results/nanoplot_after_trim/NanoPlot_CHECK_2/NanoStats.txt"

ASSEMBLY_FLYE="results/Flye_results/${species_name}_assembly.fasta"
FLYE_LOG="results/Flye_results/flye.log"
FLYE_PARAMS="results/Flye_results/params.json"

ASSEMBLY_HIFIASM="results/Hifiasm_results/${species_name}_ctg.fasta"

MAPPING_SAM="results/sam_file/${species_name}.sam"
POLISHED="results/Racon_results/${species_name}_Racon_polished.fasta"

BUSCO_FINAL=$(ls results/Busco_results/${species_name}_Busco_output/short_summary.specific.*.txt 2>/dev/null || true)
QUAST_FINAL="results/quast_report/${species_name}_Quast_output/report.txt"

############################
# 1. QC & Trimming
############################
if [[ ! -f $TRIMMED || ! -f $QC1 || ! -f $QC2 ]]; then
    echo "▶️ Running QC + trimming..."
    nextflow run 01_qc_and_trimming.nf --fastfiles "$FASTQ" --species_name "$species_name" -resume
else
    echo "QC & trimming exists, skipping."
fi

############################
# 2. Assembly (Flye or Hifiasm)
############################
if [[ "$ASSEMBLER" == "flye" ]]; then
    if [[ ! -f $ASSEMBLY_FLYE || ! -f $FLYE_LOG || ! -f $FLYE_PARAMS ]]; then
        echo "▶️ Running Flye assembly..."
        bash run_flye.sh "$TRIMMED" "$genome_size" "$flye_coverage" "$threads" "$flye_read_type" "$species_name"
    else
        echo "Flye assembly exists, skipping."
    fi
    ASSEMBLY=$ASSEMBLY_FLYE
elif [[ "$ASSEMBLER" == "hifiasm" ]]; then
    if [[ ! -f $ASSEMBLY_HIFIASM ]]; then
        echo "▶️ Running Hifiasm assembly..."
        bash run_hifiasm.sh "$species_name" "$threads"
    else
        echo "Hifiasm assembly exists, skipping."
    fi
    ASSEMBLY=$ASSEMBLY_HIFIASM
else
    echo "Unknown assembler: $ASSEMBLER"
    exit 1
fi

############################
# 3. Mapping (both assemblers)
############################
if [[ ! -f $MAPPING_SAM ]]; then
    echo "▶️ Running mapping..."
    nextflow run 03_mapping.nf \
        --fastq "$(realpath "$TRIMMED")" \
        --assembly "$(realpath "$ASSEMBLY")" \
        --species_name "$species_name" \
        -resume
else
    echo "Mapping SAM exists, skipping."
fi

############################
# 4. Polishing (Flye only)
############################
if [[ "$ASSEMBLER" == "flye" ]]; then
    if [[ ! -f $POLISHED ]]; then
        echo "▶️ Running polishing..."
        nextflow run 04_polishing.nf \
            --fastq "$TRIMMED" \
            --sam "$MAPPING_SAM" \
            --assembly "$ASSEMBLY" \
            --species_name "$species_name" \
            -resume
    else
        echo "Polishing exists, skipping."
    fi
    FINAL_ASSEMBLY=$POLISHED
else
    FINAL_ASSEMBLY=$ASSEMBLY
fi

############################
# 5. Final BUSCO + QUAST evaluation
############################
if [[ -z "$BUSCO_FINAL" || ! -f "$QUAST_FINAL" ]]; then
    echo "▶️ Running final BUSCO + QUAST evaluation..."
    nextflow run 05_eval_final.nf \
      --assembly "$FINAL_ASSEMBLY" \
      --lineage "$LINEAGE" \
      --species_name "$species_name" \
      -resume
else
    echo "Final BUSCO + QUAST evaluation exists, skipping."
fi

echo "Pipeline completed successfully."
