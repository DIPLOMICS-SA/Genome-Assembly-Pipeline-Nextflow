# Genome-Assembly-Pipeline-Nextflow
This repository contains a nextflow pipeline for denovo assembly of long nanopore reads.

# Introduction

This workflow uses the following:
* Dorado for basecalling
* Samtools for converting bam files to fastq files
* FastQC for quality check
* Nanofilt for filtering and trimming
* Flye for genome assembly
* Racon for first round assembly polishing
* Medaka for second round assembly polishing
* QUAST for assembly Quality assessment

# Dependencies

The following modules need to be loaded on the CHPC before running the pipeline:
* module purge
* module load chpc/BIOMODULES
* module load dorado
* module load samtools/1.9
* module load FastQC
* module load nanofilt
* module load flye/2.9
* module load minimap2
* module load racon/1.5.0
* module load medaka/1.11.1
* module load quast/4.6.3
* module load busco
* module load nextflow/23.10.0-all

# Usage

To obtain the workflow, having installed nextflow, users can run:
* nextflow run main.nf --help
to see the options for the workflow.

# Workflow outputs

The primary outputs of the pipeline include:
A fastq quality control report
An assembly fasta file
A quast quality report
