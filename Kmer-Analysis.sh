#!/bin/bash

######################################
#   Raw read K-Mer analysis          #
#  Wilku Meyer (wilku@cengen.co.za)  #
######################################
###################################################################################
# K-mer count with KMC and in-house R script on chpc  #
###################################################################################

file=$input
kmer=21
echo "${file}species_name_fastq_pass_con.gz" > total_number_bases.txt
zcat /mnt/lustre/users/raw/${file}species_name_fastq_pass_con.gz | awk 'NR%4==2 {sum+=length($0)} END {print sum}' >> total_number_bases.txt
kmc -cs1000 -k${kmer} -fq -t21 /mnt/lustre/users/raw/${file}species_name_fastq_pass_con.gz ${file}${kmer}mers ./                        #added -cs1000 incase there is an higher amount of coverage 
kmc_tools transform ${file}${kmer}mers histogram ${file}_kmer${kmer}_histo.txt

#############################################
# R scrict for K-mer analysis               #
# Currently runnig with Rstudio on local PC #
#############################################

R
install.packages("ggplot2")
library(ggplot2)

# Function to find peaks that are more than 10 units apart and consider 7 values for comparison and that has a higher frequency than 100 000
find_peaks <- function(data) {
  peaks <- data.frame(cov = numeric(), freq = numeric())
  for (i in 4:(nrow(data) - 3)) {
    if (data[i, "freq"] > 100000 && 
        data[i, "freq"] > data[i - 3, "freq"] &&
        data[i, "freq"] > data[i - 2, "freq"] &&
        data[i, "freq"] > data[i - 1, "freq"] &&
        data[i, "freq"] > data[i + 1, "freq"] &&
        data[i, "freq"] > data[i + 2, "freq"] &&
        data[i, "freq"] > data[i + 3, "freq"]) {
      if (nrow(peaks) == 0 || abs(data[i, "cov"] - tail(peaks$cov, 1)) > 10) {
        peaks <- rbind(peaks, data[i, ])
      }
    }
  }
  return(peaks)
}

get_ploidy <- function(peaks) {
  # Return an empty character vector if there are no peaks or the first peak's coverage is greater than 120
  if (nrow(peaks) == 0 || as.numeric(peaks$cov[1]) > 120) return(character(0))
  
  first_cov <- peaks$cov[1]  # Get the first coverage value
  
  ploidy <- sapply(peaks$cov, function(x) {
    if (x == first_cov) {
      return("n")  # Assign "n" for the first coverage value
    } else {
      ploidy_level <- round(x / first_cov)  # Calculate ploidy level
      return(paste0(ploidy_level, "n"))  # Append "n" to the calculated ploidy level
    }
  })
  
  return(ploidy)
}
# Search for the file ending with 1_histo.txt in the current directory
  file_pattern <- "_histo.txt$"
  data_file <- list.files(pattern = file_pattern)
  
  if (length(data_file) == 0) {
    warning("No file found ending with _kmer21_histo.txt in directory:", s)
    next  # Skip to the next directory
  }
  
  # Assume the total_bases file is in the same directory and has a fixed name
  total_bases_path <- 'total_number_bases.txt'
  species_info <- strsplit(data_file, "_")[[1]]
  species_name <- species_info[1]
  kmer <- species_info[2]
  
  # Read data
  data <- read.table(data_file)
  total_bases <- read.table(total_bases_path)
  colnames(data) <- c("cov", "freq")
  
  # Check if data is empty
  if (nrow(data) == 0) {
    Stop("Data is empty", s)
  }
  
  # Find peaks 
  peaks <- find_peaks(data)
  
  # Handle case when no peaks are found
  if (nrow(peaks) == 0) {
    # Create a basic plot without peaks
    y_limt <- sort(data$freq, decreasing = TRUE)
    
    p <- ggplot(data = data, aes(x = cov, y = freq)) +
      geom_line() +
      labs(x = "Coverage",
           y = "Frequency") +
      xlim(0, 250) +
      ylim(0, y_limt[2]) +  # Adjust y-axis limit
      annotate("text", x = 250, y = y_limt[2] * .95,
               label = paste(species_name, "\nNo peaks detected"," \n",kmer), hjust = 1, size = 4, fontface = "bold")
    
    # Print and save plot
    print(p)
    ggsave(filename = "kmers.png", plot = p, width = 24, height = 14, units = "cm", dpi = 300)
  } else {
    # Get ploidy levels
    ploidy_levels <- get_ploidy(peaks)
    if (length(ploidy_levels) == 0) {
      peaks$ploidy <- NA
      peaks$label <- NA
    } else {
      peaks$ploidy <- ploidy_levels
      peaks$label <- paste(peaks$ploidy, peaks$cov, sep = ":  ")
    }
    
    # Estimate genome/haploid length
    if (nrow(peaks) >1){
    length <- as.numeric(total_bases[2, 1]) / as.numeric(peaks$cov[2])
    Mb_value <- as.numeric(length) / 1e6
    formatted_value <- format(round(Mb_value, 2), big.mark = " ", decimal.mark = ".", nsmall = 2)
    length_value <- paste(formatted_value, "Mb")
    length_lable <-"Genome Length:"
    }else {length <- as.numeric(total_bases[2, 1]) / as.numeric(peaks$cov[1])
    Mb_value <- as.numeric(length) / 1e6
    formatted_value <- format(round(Mb_value, 2), big.mark = " ", decimal.mark = ".", nsmall = 2)
    length_value <- paste(formatted_value, "Mb")
    length_lable <-"Haploid Length:"
    }
    # Sort peaks by frequency in decreasing order
    y_limt <- sort(peaks$freq, decreasing = TRUE)
    
    # Adjust y_limt
    if (length(y_limt) == 0 || y_limt[1] < 1e6) {
      y_limt <- round(as.numeric(data$freq[4]))
    } else {
      y_limt <- round(as.numeric(y_limt[1]) * 1.1)
    }
    # Set x_limt based on the last peak's coverage
    last_peak <- sort(peaks$cov, decreasing = TRUE)
    if (length(last_peak) > 0 && 250 < as.numeric(last_peak[1]) && as.numeric(last_peak[1]) < 1000) {
      x_limt <- as.numeric(last_peak[1]) + 100
    } else {
      x_limt <- 250
    }
    
    # Create ggplot
    p <- ggplot(data = data, aes(x = cov, y = freq)) +
      geom_line() +
      labs(x = "Coverage",
           y = "Frequency") +
      xlim(0, x_limt) +
      ylim(0, y_limt) +
      geom_vline(data = peaks, aes(xintercept = cov), linetype="dotted") +
      geom_text(data = peaks, aes(x = cov, y = freq, label = label),
                vjust = -0.5, hjust = 0.5, size = 4) +
      annotate("text", x = x_limt, y = y_limt * 0.95,
               label = paste(species_name, "\nEstimated ",length_lable, length_value,"\n ",kmer), hjust = 1, size = 4, fontface = "bold")
    p
    ggsave(filename = "kmers.png", plot = p, width = 24, height = 14, units = "cm", dpi = 300)
  }
  
