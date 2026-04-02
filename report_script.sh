#!/bin/bash
set -euo pipefail
########################

#Usage: bash report_script.sh species_name 

#species_name should be the same as in the params.config file

########################

# === USER INPUT ===
species_name="${1:-}"

if [[ -z "$species_name" ]]; then
    echo "❌ Usage: bash report_script.sh <species_name>"
    exit 1
fi

cd ./results

# === Load dependencies ===
module load chpc/BIOMODULES
module load samtools

################################################################################
# PART 1 — SAMTOOLS PROCESSING
################################################################################
echo "=== PART 1: SAMTOOLS ==="

if [[ ! -f ./sam_file/${species_name}.sam ]]; then
    echo "❌ ${species_name}.sam not found in ./sam_file/"
    exit 1
fi

if [[ ! -f ./sam_file/${species_name}_sorted.bam ]]; then
    echo "▶️ Sorting SAM..."
    samtools view -bS ./sam_file/${species_name}.sam | \
        samtools sort -o ./sam_file/${species_name}_sorted.bam
else
    echo "✅ Sorted BAM already exists"
fi

if [[ ! -f ./sam_file/minimap2_coverage.txt ]]; then
    echo "▶️ Calculating depth coverage..."
    samtools depth ./sam_file/${species_name}_sorted.bam |
        awk '{sum+=$3} END { print "Average = ",sum/NR}' \
        > ./sam_file/minimap2_coverage.txt
else
    echo "✅ Coverage already calculated"
fi

if [[ ! -f ./sam_file/sam_stats.txt ]]; then
    echo "▶️ Generating SAM stats..."
    samtools stats ./sam_file/${species_name}_sorted.bam |
        grep ^SN | cut -f 2- > ./sam_file/sam_stats.txt
else
    echo "✅ SAM stats already exist"
fi

################################################################################
# PART 2 — Detect Assembler
################################################################################
echo "=== PART 2: Assembler Detection ==="

assembler=""
if [[ -d ./Flye_results ]]; then
    assembler="flye"
elif [[ -d ./Hifiasm_results ]]; then
    assembler="hifiasm"
fi

if [[ "$assembler" == "flye" ]]; then
    echo "🔎 Flye assembly detected — extracting mean coverage..."
    OUTPUT_FILE="${species_name}_flye_mean_coverage.txt"
    > "$OUTPUT_FILE"
    for logfile in $(find ./Flye_results -type f -name "flye.log"); do
        mean_cov=$(grep "Mean coverage:" "$logfile" | awk '{print $3}')
        [[ -n "$mean_cov" ]] && echo "$logfile: $mean_cov" >> "$OUTPUT_FILE"
    done
elif [[ "$assembler" == "hifiasm" ]]; then
    echo "🔎 Hifiasm assembly detected — no flye.log coverage to extract."
else
    echo "⚠️ No assembler results directory found."
fi

################################################################################
# PART 3 — Collect Outputs
################################################################################
echo "=== PART 3: Collect Outputs ==="

mkdir -p "${species_name}" "${species_name}_other_results_outputs"

mv ./nanoplot_before_trim/NanoPlot_CHECK_1/NanoStats.txt "${species_name}_NanoStats_before_trim.txt" 2>/dev/null || true
mv ./nanoplot_before_trim/NanoPlot_CHECK_1/NanoPlot-report.html "${species_name}_NanoPlot_before_trim.html" 2>/dev/null || true
mv ./nanoplot_after_trim/NanoPlot_CHECK_2/NanoStats.txt "${species_name}_NanoStats_after_trim.txt" 2>/dev/null || true
mv ./nanoplot_after_trim/NanoPlot_CHECK_2/NanoPlot-report.html "${species_name}_NanoPlot_after_trim.html" 2>/dev/null || true

mv ./trimmed_fastq/${species_name}.trimmed.fastq ./"${species_name}" 2>/dev/null || true

if [[ "$assembler" == "flye" ]]; then
    mv ./Flye_results/${species_name}_assembly.fasta ./"${species_name}/${species_name}_flye_assembly.fasta" 2>/dev/null || true
    mv ./Racon_results/${species_name}_Racon_polished.fasta ./"${species_name}/${species_name}_racon_polished.fasta" 2>/dev/null || true
elif [[ "$assembler" == "hifiasm" ]]; then
    mv ./Hifiasm_results/${species_name}_ctg.fasta ./"${species_name}/${species_name}_hifiasm_assembly.fasta" 2>/dev/null || true
fi

mv ./Busco_results/Busco_output/short_summary.specific.*.txt "${species_name}_busco_summary.txt" 2>/dev/null || true
mv ./quast_report/Quast_result/report.txt "${species_name}_quast_report.txt" 2>/dev/null || true

mv ./sam_file/minimap2_coverage.txt "${species_name}_minimap2_coverage.txt"
mv ./sam_file/sam_stats.txt "${species_name}_sam_stats.txt"

################################################################################
# PART 4 — Final Report
################################################################################
echo "=== PART 4: Final Report ==="

ordered_files=(
    "${species_name}_kmer_total_number_bases.txt"
    "${species_name}_kmer_cov_size.txt"
    "${species_name}_NanoStats_before_trim.txt"
    "${species_name}_NanoStats_after_trim.txt"
    "${species_name}_minimap2_coverage.txt"
    "${species_name}_sam_stats.txt"
    "${species_name}_quast_report.txt"
    "${species_name}_busco_summary.txt"
)

[[ -f "${species_name}_flye_mean_coverage.txt" ]] && ordered_files+=("${species_name}_flye_mean_coverage.txt")

report="${species_name}_report.txt"
{
    echo "${species_name^} Assembly Report for RedCap"
    echo "Generated on: $(date)"
    echo -e "\n========================================\n"
    for file in "${ordered_files[@]}"; do
        [[ -f "$file" ]] && { echo "========== $file =========="; cat "$file"; echo -e "\n\n"; }
    done
} > "$report"

echo "✅ Report generated: $report"

################################################################################
# PART 5 — Organize Outputs
################################################################################
echo "=== PART 5: Organizing Outputs ==="

mv *_NanoStats_before_trim.txt *_NanoStats_after_trim.txt *_NanoPlot_*.html "${species_name}" 2>/dev/null || true
mv "${species_name}"_busco_summary.txt "${species_name}" 2>/dev/null || true
mv "${species_name}"_quast_report.txt "${species_name}" 2>/dev/null || true

mv "${species_name}_minimap2_coverage.txt" "${species_name}_other_results_outputs" 2>/dev/null || true
mv "${species_name}_sam_stats.txt" "${species_name}_other_results_outputs" 2>/dev/null || true
mv "${species_name}_flye_mean_coverage.txt" "${species_name}_other_results_outputs" 2>/dev/null || true

mv nanoplot_before_trim nanoplot_after_trim quast_report Busco_results trimmed_fastq sam_file \
   Flye_results Racon_results Hifiasm_results "${species_name}_other_results_outputs" 2>/dev/null || true

################################################################################
# PART 6 — Software Versions
################################################################################
echo "=== PART 6: Appending Software Versions ==="

CONFIG="../params.config"
[[ ! -f "$CONFIG" ]] && { echo "❌ params.config not found"; exit 1; }

source "$CONFIG"
species_name="$1"  # ensure CLI takes priority

module load chpc/singularity

CONTAINER="/home/apps/chpc/bio/1ksa_pipeline/1ksa_pipeline.sif"
BUSCO_CONTAINER="/home/apps/chpc/bio/busco/5.8.0/busco_v5.8.0_cv1.sif"

CONTAINER_TOOLS=(NanoPlot NanoFilt hifiasm flye minimap2 samtools quast)
DIRECT_MODULE_TOOLS=(chopper "kraken2/2.1.2" "racon/1.5.0")

source /etc/profile.d/modules.sh
module load chpc/BIOMODULES

{
    echo -e "\n========================================"
    echo "Software Versions"
    echo "Checked on: $(date)"
    echo "========================================"
} >> "$report"

# Container tools
for TOOL in "${CONTAINER_TOOLS[@]}"; do
    echo ">> $TOOL version:" >> "$report"
    singularity exec "$CONTAINER" $TOOL --version 2>&1 \
        | grep -v "WARNING" | grep -v "^$" >> "$report" || true
    echo "" >> "$report"
done

# BUSCO container
echo ">> BUSCO version:" >> "$report"
singularity exec "$BUSCO_CONTAINER" busco --version 2>&1 \
    | grep -v "WARNING" | grep -v "^$" >> "$report" || true

# Module tools
for TOOL in "${DIRECT_MODULE_TOOLS[@]}"; do
    DISPLAY="${TOOL%%/*}"
    echo ">> $DISPLAY version:" >> "$report"

    if module load "$TOOL" 2>/dev/null; then
        if [[ "$DISPLAY" == "racon" ]]; then
            echo "   racon version: ${TOOL##*/}" >> "$report"
        else
            $DISPLAY --version 2>&1 | head -1 >> "$report" || true
        fi
    fi
    echo "" >> "$report"
done

echo "✅ Software versions appended"

echo "🎉 Done!"
echo "Report: $report"
