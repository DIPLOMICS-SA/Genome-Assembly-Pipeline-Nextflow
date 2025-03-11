# 1KSA-Genome-Assembly-Pipeline-Nextflow
## Project Overview
This guide provides instructions for using the 1KSA Nextflow pipeline to perform de novo genome assembly with Oxford Nanopore long-read sequencing data. Please note that the de novo assembly process generates a draft-level assembly.
Purpose: This project supports the 1KSA initiative, focusing on sequencing and assembling indigenous South African species. The draft-level assembly serves as an essential layer of quality control, demonstrating that there are sufficient, high-quality data that can be used as a tool for biodiversity genomics. By contributing to biodiversity research, this work aims to tackle critical global challenges such as reducing food scarcity and protecting native species. Advances in biodiversity genomics empower researchers to combat threats to biodiversity and maintain genetic diversity. The initial phase of biodiversity research involves identifying and analyzing species to better understand ecosystems, their functions, and interdependencies, highlighting the importance of building reference genomes.
1KSA website: https://www.1ksa.org.za
Detailed instructions: https://zenodo.org/communities/1ksa/records?q=&l=list&p=1&s=10&sort=newest

# Workflow Components

The 1KSA Genome Assembly Pipeline (built in nextflow) is intended for use on the Centre for High Performance Computing (CHPC) and uses the following tools (Figure 1):
* KMC for counting of k-mers in DNA (done separately)
* Samtools for converting bam files to fastq files
* Nanoplot for quality check
* Nanofilt for filtering and trimming
* Flye v2.9 for genome assembly
* Racon v1.5.0 for assembly polishing
* BUSCO v5.4.5 for assembly quality assessment
* QUAST v5.2.0 for assembly quality assessment

The starting point for this workflow is RAW fastq files, i.e. basecalling has already been done.

![Image Alt text](https://github.com/DIPLOMICS-SA/Genome-Assembly-Pipeline-Nextflow/blob/f8a896f93a1469db564ac9fc4a78bd26db588d3b/1KSA_assmbly_pipeline_07-03-2025.png)
Figure 1: 1KSA Workflow 


# 1. CHPC login and Pipeline Download
1.1 Login to the CHPC using your own user account

Link to CHPC quick start guide: https://wiki.chpc.ac.za/quick:start

``` ssh username@lengau.chpc.ac.za  ```


The following modules need to be loaded on the CHPC before running the pipeline:

```
module purge
module load chpc/BIOMODULES
module load dorado
module load samtools/1.9
module load nanoplot
module load nanofilt
module load flye/2.9
module load minimap2
module load racon/1.5.0
module load medaka/1.11.3
module load quast/4.6.3
module load busco/5.4.5
module load bbmap/38.95
module load metaeuk
module load python
module load R
module load KMC
module load genomescope
module load smudgeplot
module load nextflow/23.10.0-all
```
The following models and databases need to be downloaded before running the pipeline:
* Dorado:

``` dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0 ```

* Busco:

``` busco --download eukaryota_odb10 ```

# Usage

To obtain the workflow, having installed nextflow, users can run:
* nextflow run main.nf --help
to see the options for the workflow.

# Workflow outputs

The primary outputs of the pipeline include:
* A fastq quality control report
* A busco report before polishing 
* A quast quality report before polishing 
* 2 assembled fasta files (From Flye and Racon)
* A busco report after polishing - The BUSCO score from this report is documented on the species card
* A quast quality report after polishing - The genome size from this report is documented on the species card
