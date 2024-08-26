###############################################################
#   Raw read K-Mer analysis                                   #
#  Wilku Meyer (wilku@cengen.co.za)                           #
#  Setshaba Taukobong (setshaba.taukobong@diplomics.org.za)   #
###############################################################
###################################################################################
# K-mer count with KMC on chpc, Jellyfish will also work but KMC is much quicker  #
###################################################################################

#!/bin/bash

# Add KMC module
module load KMC

set -e

file=Species_name.fastq
kmer=21

echo "${file}" > total_number_bases.txt
zcat ${file} | awk 'NR%4==2 {sum+=length($0)} END {print sum}' >> total_number_bases.txt
kmc -k${kmer} -fq -t20 ${file} ${kmer}mers .
kmc_tools transform ${kmer}mers histogram kmer${kmer}.histo

# Add genomescope-2 module
module load genomescope-2

genomescope.R -i kmer${kmer}.histo -o GenomeScope-Results -k${kmer}


### Double check the ploidy level of your species
# Add smudgeplot module

L=$(smudgeplot.py cutoff kmer${kmer}.histo L)
U=$(smudgeplot.py cutoff kmer${kmer}.histo U)
echo $L $U 
kmc_tools transform ${kmer}mers -ci"$L" -cx"$U" dump -s ${kmer}mers_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o ${kmer}mers_L"$L"_U"$U" < ${kmer}mers_L"$L"_U"$U".dump
smudgeplot.py plot ${kmer}mers_L"$L"_U"$U"_coverages.tsv
