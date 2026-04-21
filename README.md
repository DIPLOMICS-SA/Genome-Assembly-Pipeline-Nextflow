# 1KSA-Genome-Assembly-Pipeline-Nextflow
## Project Overview
This guide provides instructions for using the 1KSA Nextflow pipeline to perform de novo genome assembly with Oxford Nanopore long-read sequencing data. Please note that the de novo assembly process generates a draft-level assembly.

Purpose: This project supports the 1KSA initiative, focusing on sequencing and assembling indigenous South African species. The draft-level assembly serves as an essential layer of quality control, demonstrating that there are sufficient, high-quality data that can be used as a tool for biodiversity genomics. By contributing to biodiversity research, this work aims to tackle critical global challenges such as reducing food scarcity and protecting native species. Advances in biodiversity genomics empower researchers to combat threats to biodiversity and maintain genetic diversity. The initial phase of biodiversity research involves identifying and analyzing species to better understand ecosystems, their functions, and interdependencies, highlighting the importance of building reference genomes.

1KSA website: https://www.1ksa.org.za

Detailed instructions: [https://zenodo.org/communities/1ksa/records?q=&l=list&p=1&s=10&sort=newest](https://doi.org/10.5281/zenodo.17811930)

## Decision Tree
Once sequencing is complete and your data has been transferred to the CHPC, refer to the decision tree below to determine the most suitable assembly approach moving forward:
![Image Alt text](https://github.com/DIPLOMICS-SA/Genome-Assembly-Pipeline-Nextflow/blob/main/Figure_1.png)
Figure 1: 1KSA Draft genome assembly decision tree

## Workflow Components:

The 1KSA Genome Assembly Pipeline is intended for use on the Centre for High Performance Computing (CHPC) and uses the following tools (Figure 2):
* KMC for counting of k-mers in DNA (done separately)
* Nanoplot v1.46.1 for quality check
* Chopper v0.10.0 for filtering and trimming
* Flye v2.9.6 for genome assembly (genomes smaller than 3 Gb)
* Racon v1.5.0 for assembly polishing (flye assembly only)
* Hifiasm v0.25.0 for genome assembly (genomes larger than 3 Gb)
* BUSCO v6.0.0 for assembly quality assessment
* QUAST v5.3.0 for assembly quality assessment

The starting point for this workflow is RAW fastq files, i.e. basecalling has already been done.

![Image Alt text](https://github.com/DIPLOMICS-SA/Genome-Assembly-Pipeline-Nextflow/blob/9c2401b6741398cc1cb165bd63f13da544b7b47e/Figure_2.png)

Figure 2: 1KSA Workflow


## 1. CHPC login and Pipeline Download
1.1 Login to the CHPC using your own user account

Link to CHPC quick start guide: https://wiki.chpc.ac.za/quick:start

``` 
ssh username@lengau.chpc.ac.za
```

1.2 Clone the Pipeline Repository with all the necessary scripts:

* kmer-Analysis.sh
* kmerPlot.R
* master.sh
* params.config
* run_flye.sh
* run_hifiasm.sh
* nextflow.config
* 01_qc_and_trimming.nf
* 03_mapping.nf
* 04_polishing.nf
* 05_eval_final.nf

```
## Navigate to the directory with your fastq data
cd /path/to/folder/with/species/fastq/files                       
              
git clone https://github.com/DIPLOMICS-SA/Genome-Assembly-Pipeline-Nextflow.git

## Navigate inside the folder with the pipeline:
cd Genome-Assembly-Pipeline-Nextflow
```

1.3 Download the BUSCO Lineage Folder
```
module load chpc/BIOMODULES
module load busco/5.8.0
singularity run $SIF busco --download eukaryota_odb10             #change database if necessary (viridiplantae_odb10; insecta_odb10)

## Make sure the download was completed:
ls ./busco_downloads/lineages/eukaryota_odb10                     #change database if necessary (viridiplantae_odb10; insecta_odb10)
```

## 2. CHPC PBS Session setup

For this pipeline, resources will be requested from two queues on the CHPC: seriallong for initial quality control and bigmem for assembly and polishing.

* Seriallong: 128 GB RAM, max wall time is 144 hours
* Bigmem: 1 TB RAM, max wall time is 48 hours

This workflow requires 3 screen sessions (Figure 2):

* screen_1 (seriallong, 100 hours)
    * Unzip and concatenate fastq.gz files
    *  Perform quality control (kmer analysis and filtering)

* screen_2 (bigmem, 48 hours)
    * Genome assembly and polishing using Flye and Racon, respectively

* screen_3 (seriallong, 100 hours)
    * Evaluation with Quast and BUSCO

![Image Alt text](https://github.com/DIPLOMICS-SA/Genome-Assembly-Pipeline-Nextflow/blob/169aec45b0de1218a81849ab5a806afc4d053029/Figure_2.png)

Figure 3: 1KSA Workflow processes in different CHPC queues depending on capacity requirements

## 3. Data Preparation
### 3.1 Screen 1 (Quality Control)
3.1.1 Start an interactive job in screen 1

```
screen -S screen_1

qsub -I -l select=1:ncpus=12:mpiprocs=1 -q seriallong -P CBBIXXXX -l walltime=100:00:00   #change project number
```

3.1.2 Concatenate FASTQ Files

```
## If you have not concatenated the fastq files, navigate to the folder (pass) and concatenate it:
cd /path/to/folder/with/species/fastq/files                       #change the path
cat *.fastq > species_name_fastq_pass_con.fastq                   #change the species_name
```

3.1.3 Unzip Files (If Required)

```
## If you need to unzip the files first, do the following:
gunzip -d *.gz | cat *.fastq > species_name_fastq_pass_con.fastq  #change the species_name
```

```
## If you struggle with file permissions:
cat ./fastq_pass/*.gz > species_name_fastq_pass_con.fastq.gz
gunzip -d species_name_fastq_pass_con.fastq.gz
```

## 4. K-mer Analysis

### 4.1 Screen 1 (Quality Control)

This step is crucial for determining the read coverage and estimating the genome size, both of which are essential parameters for optimizing the Flye assembly process. Two scripts are used for this step, which is included in the Genome-Assembly-Pipeline-Nextflow repository:

* kmer-analysis.sh
* kmerPlot.R (used by the first script)

This first script loads the necessary modules, calculates k-mer frequencies from a FASTQ file using KMC, generates a k-mer histogram, and runs an R script for k-mer frequency analysis on raw sequencing reads to estimate genome size, coverage, and ploidy.

#### What it does:

##### 1. K-mer Counting:
* Uses KMC to count 21-mers (or other specified k-mer sizes) from input fastq reads.
* Outputs a histogram of k-mer coverage vs frequency.

##### 2. Total Bases Calculation:
* Calculates the total number of bases from the input fastq file (sum of all sequence lengths).

##### 3. Peak Detection:
* Identifies coverage peaks in the k-mer histogram, which correspond to different ploidy levels (haploid, diploid, etc.).
* Detects plateaus in the histogram to account for repetitive regions.

##### 4. Genome Size Estimation:
* Estimates genome size based on the position of peaks and total bases.
* Outputs genome size in megabases (Mb).

##### 5. Ploidy Inference:
* Labels peaks as n (haploid), 2n (diploid), etc., based on their relative coverage.

##### 6. Visualization:
* Generates plots (PDF and PNG) of the k-mer histogram with annotated peaks and plateaus.

##### Output files:
* total_number_bases.txt
* species_name_17_mers.kmc_suf
* species_name_17_mers.kmc_pre
* species_name_17_mers_histo.txt
* species_name_k_mers.png
* species_name_k_mers.pdf

#### 4.1.1 Run K-mer Analysis (continue in the interactive session in screen_1)

```
## Navigate to your working directory
cd /path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow
#change path 

bash kmer-Analysis.sh /path/to/fastq/file/species_name_fastq_pass_con.fastq species_name
#change path
#change species_name
```

#### 4.1.2 Record Genome Coverage and Size

```
zless k_mers_Stats_species_name.txt
```
Kmer analysis example output:
```
Species_name: k-mer= 17
Total input bases 153410333621
Peaks or Plateaus detected=2
Ploidy= 2n =2n
2n Genome Length=0.87 Gb at 176 X Coverage
Expected Assembly Length if fully collapsed=0.87 Gb at 176 X Coverage
```

You will use this as input parameters when running the pipeline below.

Download the .png and view on your local computer.
```
# Open a new terminal and navigate to the Desktop
cd Desktop
# Copy the .png file to your local computer using rsync
rsync -av username@lengau.chpc.ac.za:/path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow/species_name_k_mers.png .
# change username
# Enter password
```

## 5. Run the Nextflow Pipeline

### 5.1 Screen 1 (Quality Control)

#### 5.1.1 Re-attach to screen 1

```
screen -r screen_1

## Navigate to your working directory
cd /path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow
#change path
```

#### 5.1.2 Open the params.config file:

```
nano params.config
```

#### 5.1.3: Edit the params.config file

Change the following values:

```
# File: params.config

species_name='species_name'    # Change species_name
assembler='hifiasm'            # Options: 'flye' or 'hifiasm'

# Shared parameters
threads=15
LINEAGE='eukaryota_odb10'      #Mandatory BUSCO lineage (Other options: insecta_odb10, viridiplantae_odb10)

# Flye-specific parameters (mandatory)
genome_size='0.77g'            # Mandatory genome size - adjust to the values obtained from the k-mer analysis
flye_coverage=45               # Mandatory coverage - adjust to the values obtained from the k-mer analysis
flye_read_type='nano-raw'      # Default read type (can change to nano-hq for Q15 reads)
```

#### 5.1.4 Save the script
```
## Save and Exit
^X  # Exit  
y   # Confirm Save  
Enter
```

#### 5.1.5 Run the pipeline:

```
module load chpc/BIOMODULES
module load nextflow/24.04.4-all
bash master.sh /path/to/fastq/file      # Change path
```

#### Guidelines for Adjusting Assembly Parameters
When running the pipeline, modify the following parameters based on your k-mer analysis and sequencing read quality:

##### 1. Coverage and Genome Size
* Set ```flye_coverage``` and ```genome_size``` according to the k-mer analysis output.

##### 2. Read Type Selection
* Choose nano-raw for low-quality (Q > 10) nanopore reads.
* Choose nano-hq for high-quality (Q > 15) nanopore reads.

##### 3. BUSCO Lineage Selection
* Select an appropriate BUSCO lineage dataset based on your organism:
* Use ```eukaryota_odb10``` for general eukaryotic genomes.
* Use a more specific dataset based on your organism type, such as, ```viridiplantae_odb10``` for plants, or insecta_odb10 for insects.

#### Table 1: Genome assembly pipeline parameters
| Tool     | Parameter Description          | Parameter Flag       | Value/Notes                 |
|----------|--------------------------------|----------------------|-----------------------------|
| NanoFilt | Quality                        | -q                   | 10                          |
|          | Minimum read length            | -l                   | 1000                        |
| Flye     | ONT regular reads (Q>10)       | --nano-raw           |                             |
|          | ONT high-quality reads (Q>15)  | --nano-hq            |                             |
|          | Threads                        | -t                   | 15                          |
|          | Reduced coverage               | --asm-coverage       | (k-mer analysis) e.g. 36    |
|          | Genome size                    | -g                   | (k-mer analysis) e.g. 2.4g  |
| Racon    | Match                          | -m                   | 3                           |
|          | Mismatch                       | -x                   | -5                          |
|          | Gap                            | -g                   | -4                          |
|          | Window length                  | -w                   | 500                         |
|          | Threads                        | -t                   | 15                          |
| BUSCO    | Database                       | -l                   | eukaryota_odb10             |
|          | Mode                           | -m                   | genome                      |
|          | Metaeuk gene predictor         | --metaeuk_parameters | METAEUK_PARAMETERS --offline|


Detach from screen_1: ```CRTL A+D```

### 5.2 Screen 2 (Assembly & Polishing)

When the pipeline runs out of memory (at the Flye step), switch over to a bigmem node and resume the pipeline using the exact same command you used to start the pipeline: bash master.sh /path/to/fastq/file (screen_2).

#### 5.2.1 Start the second screen session

```
screen -S screen_2

## Navigate to your working directory
cd /path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow      #change path
```

#### 5.2.2 Create a Job Script (bigmem)

```
nano name_of_your_scrip_1.sh                          #change the name of the script
```

#### 5.2.3 Example Script Template (copy and paste)

```
#!/bin/bash
#PBS -l select=1:ncpus=40                             #change ncpus if necessary
#PBS -l walltime=48:00:00
#PBS -q bigmem
#PBS -P CBBIXXXX                                      #change project number
#PBS -o /path/to/your/working/directory/stdout.txt    #change path
#PBS -e /path/to/your/working/directory/stderr.txt    #change path
#PBS -M email@address.org.za                          #change email
#PBS -m b

#############################################################################################
cd /path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow   # change path

fastq='/path/to/fastq/file'                                                     # change path
#############################################################################################

module load chpc/BIOMODULES
module load nextflow/24.04.4-all

bash master.sh $fastq

```

#### 5.2.4 Save and submit the script

```
## Save and Exit
^X  # Exit  
y   # Confirm Save  
Enter  

## Submit job script  
qsub name_of_your_script_1.sh
```

#### 5.2.5 If the wall time runs out, resubmit your scrip:

```
## Navigate to your working directory
cd /path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow    #change path

## Submit your job script:  
qsub name_of_your_scrip_1.sh
```

Detach from screen_2: ```CRTL A+D```

### 5.3 Screen 3 (Evaluation)

#### 5.3.1 If the wall time runs out and the assembly is complete, start the screen 3 session

```
screen -S screen_3

## Navigate to your working directory
cd /path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow  
```

#### 5.3.2 Create a job script (seriallong)

```
nano name_of_your_scrip_2.sh
```

#### 5.3.3 Example Script Template (copy and paste)

```
#!/bin/bash
#PBS -l select=1:ncpus=12                             #change ncpus if necessary
#PBS -l walltime=100:00:00
#PBS -q seriallong
#PBS -P CBBIXXXX                                      #change project number
#PBS -o /path/to/your/working/directory/stdout.txt    #change path
#PBS -e /path/to/your/working/directory/stderr.txt    #change path
#PBS -M email@address.org.za                          #change email
#PBS -m b

#############################################################################################
cd /path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow   # change path

fastq='/path/to/fastq/file'                                                     # change path
#############################################################################################

module load chpc/BIOMODULES
module load nextflow/24.04.4-all

bash master.sh $fastq
```

#### 5.3.4 Save and submit the script

```
## Save and Exit
^X  # Exit  
y   # Confirm Save  
Enter  

## Submit job script  
qsub name_of_your_script_2.sh
```

Detach from screen_3: ```CRTL A+D```

## Generate a report of the assembly metrics
### 6.1 Screen 4

#### 6.1.1 Start an interactive job in screen 4

```
##  Re-attach to screen 4
screen -r screen_4

## Start an interactive job:
qsub -I -l select=1:ncpus=12:mpiprocs=1 -q seriallong -P CBBIXXXX -l walltime=100:00:00   #change project number

## Navigate to your working directory
cd /path/to/folder/with/species/fastq/files/Genome-Assembly-Pipeline-Nextflow
```

#### 6.1.2 Rename the assembly output files and generate a single report

```
## Run the report script:
bash report_script.sh species_name
#change the species_name
```

# Workflow outputs

The primary outputs of the pipeline include:
* A fastq quality control report
* A busco report before polishing 
* A quast quality report before polishing 
* 2 assembled fasta files (from Flye and Racon)/ 1 assembled fasta file (from Hifiasm)
* A busco report after polishing - The BUSCO score from this report is documented on the species card
* A quast quality report after polishing - The genome size from this report is documented on the species card
