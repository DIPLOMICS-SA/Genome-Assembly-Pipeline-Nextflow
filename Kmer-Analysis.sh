######################################
#   Raw read K-Mer analysis          #
#  Wilku Meyer (wilku@cengen.co.za)  #
######################################
###################################################################################
# K-mer count with KMC on chpc, Jellyfish will also work but KMC is much quicker  #
###################################################################################

#!/bin/bash

set -e

file=Spekboom_unfragmented.fastq
kmer=21

echo "${file}" > total_number_bases.txt
cat "${file}" | awk 'NR%4==2 {sum+=length($0)} END {print sum}' >> total_number_bases.txt
kmc -k${kmer} -fq -t21 ${file} ${kmer}mers .
kmc_tools transform ${kmer}mers histogram kmer${kmer}_histo.txt

# Add R module (includes appropriate openMPI and gcc modules)
module load R/4.2.0

# Run R script
Rscript kmerPlot.R
