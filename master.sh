#!/bin/bash
set -euo pipefail

# Load parameters
CONFIG="params.config"
if [[ ! -f $CONFIG ]]; then
    echo "❌ Missing $CONFIG file. Please create it first."
    exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG"

# Positional argument: raw FASTQ file
FASTQ=${1:?Please provide FASTQ input file}

# Load lineage from config (already sourced)
LINEAGE="$lineage"

# Derived paths
# ---------- Trimmed output ----------
TRIMMED="results/trimmed_fastq/sample_id.trimmed.fastq"
QC1="results/nanoplot_before_trim/NanoPlot_CHECK_1/NanoStats.txt"
QC2="results/nanoplot_after_trim/NanoPlot_CHECK_2/NanoStats.txt"

# ---------- Flye assembly output ----------
ASSEMBLY="results/Flye_results/assembly.fasta"
FLYE_LOG="results/Flye_results/flye.log"
FLYE_PARAMS="results/Flye_results/params.json"

# ---------- BUSCO before polishing ----------
BUSCO1=$(ls results/Busco_results/Busco_outputs1/short_summary.specific.*.txt 2>/dev/null || true)

# ---------- QUAST before polishing ----------
QUAST1="results/quast_report/Quast_output1/report.txt"

# ---------- Polishing outputs ----------
POLISHED="results/Racon_results/Racon_polished.fasta"

# ---------- BUSCO + QUAST after polishing ----------
BUSCO_POLISHED=$(ls results/Busco_results/Busco_outputs2/short_summary.specific.*.txt 2>/dev/null || true)
QUAST2="results/quast_report/Quast_output2/report.txt"

############################
# 1. QC & Trimming
############################
if [[ ! -f $TRIMMED || ! -f $QC1 || ! -f $QC2 ]]; then
    echo "▶️  Running QC + trimming..."
    nextflow run 01_qc_and_trimming.nf --fastfiles "$FASTQ" -resume
else
    echo "✅  QC & trimmiing complete, skipping."
fi

############################
# 2. Flye Assembly
############################
if [[ ! -f $ASSEMBLY || ! -f $FLYE_LOG || ! -f $FLYE_PARAMS ]]; then
    echo "▶️  Running Flye..."
    bash run_flye.sh "$TRIMMED" "$genome_size" "$flye_coverage" "$flye_threads" "$flye_read_type"
else
    echo "✅  Flye assembly exists, skipping."
fi

############################
# 3. Polishing & Evaluation
############################
if [[ -z "$BUSCO1" || ! -f "$QUAST1" ]]; then
    echo "▶️  Running assembly evaluation (pre-polish)..."
    nextflow run 01_eval_before_polishing.nf \
        --assembly "$ASSEMBLY" \
        --lineage "$lineage" \
        -resume
else
    echo "✅  Pre-polishing evaluation (BUSCO + QUAST) complete, skipping."
fi

if [[ ! -f $POLISHED || -z "$BUSCO_POLISHED" || ! -f $QUAST2 ]]; then
    echo "▶️  Running polishing and evaluation..."
    exec nextflow run 02_polishing_and_eval.nf \
      --fastfiles "$FASTQ" \
      --assembly "$ASSEMBLY" \
      --lineage "$LINEAGE" \
      -resume
else
    echo "✅  Polishing complete, skipping."
fi

echo "✅ Pipeline completed successfully."

