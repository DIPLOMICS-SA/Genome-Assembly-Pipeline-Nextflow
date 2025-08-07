######################################
#   Raw read K-Mer analysis          #
#  Wilku Meyer (wilku@cengen.co.za)  #
######################################
###################################################################################
# K-mer count with KMC on chpc, Jellyfish will also work but KMC is much quicker  #
###################################################################################

#!/bin/bash

set -e

# Add R module (includes appropriate openMPI and gcc modules)
module load R
module load KMC

# Check if input file is provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_file> <genus_species>"
    exit 1
fi

# Input file
file=$1
kmer=17  # Set the k-mer size here
species=$2 
echo "${file}" > total_number_bases.txt
cat "${file}" | awk 'NR%4==2 {sum+=length($0)} END {print sum}' >> total_number_bases.txt &
kmc  -cs1000 -m64 -sm -ci3 -k${kmer} -fq -t7 ${file} ${species}_${kmer}_mers .
kmc_tools transform ${species}_${kmer}_mers histogram ${species}_${kmer}_mers_histo.txt
wait
# Run R script
Rscript kmerPlot.R