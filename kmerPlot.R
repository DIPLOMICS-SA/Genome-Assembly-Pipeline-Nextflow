##################################################
#R script for K-mer analysis
#Wilku Meyer 
#contact: wilku@cengen.ac.za
##################################################

# Function to find peaks that are more than 10 units apart and consider 7 values for comparison
find_peaks <- function(data) {
  peaks <- data.frame(cov = numeric(), freq = numeric())
  for (i in 4:(nrow(data) - 3)) {
    if (data[i, "cov"] < 1000 &&
        data[i, "freq"] > 500000 &&
        data[i, "freq"] > data[i - 3, "freq"] &&
        data[i, "freq"] > data[i - 2, "freq"] &&
        data[i, "freq"] > data[i - 1, "freq"] &&
        data[i, "freq"] > data[i + 1, "freq"] &&
        data[i, "freq"] > data[i + 2, "freq"] &&
        data[i, "freq"] > data[i + 3, "freq"]  
    ) {
      if (nrow(peaks) == 0 || abs(data[i, "cov"] - tail(peaks$cov, 1)) > 10) {
        peaks <- rbind(peaks, data[i, ])
      }
    }
  }
  return(peaks)
}

# Function to find plateaus that are more than 10 units apart and are 3 values or longer
find_plateaus <- function(data) {
  plateaus <- data.frame(cov = numeric(), freq = numeric())
  plateaus_end <- data.frame(cov = numeric(), freq = numeric())
  n <- nrow(data)
  i <- 1
  while (i <= (n - 1)) {
    start <- i
    end <- i
    diffs <- numeric()
    
    # Collect differences while within the range and values are positive
    while (end < n &&
           data[end, "freq"] > 1000000) {
      
      diff <- (data[end , "freq"] - data[end+ 1, "freq"]) / data[end, "freq"]
      
      if (diff > -0.0075 &&diff < 0.0075) { #  Small difference threshold for plateaus
        diffs <- c(diffs, diff)
        end <- end + 1
      } else {
        break
      }
    }
    
    # Check conditions for plateau
    if (length(diffs) >= 3 &&
        (nrow(plateaus) == 0 || abs(data[start, "cov"] - tail(plateaus_end$cov, 1)) > 10)
    ) { 
      
      # Store the start of the plateau
      plateaus <- rbind(plateaus, data[round((start+end)/2), , drop = FALSE])
      plateaus_end <- rbind(plateaus, data[end, , drop = FALSE])
    }
    
    # Move to the next potential plateau
    i <- end + 1
  }
  
  return(plateaus)
}

# Custom merge function
merge_peaks_and_plateaus <- function(peaks, plateaus) {
  merged <- peaks
  
  # Determine the threshold
  if(is.na(peaks$cov[1]) || peaks$cov[1] < 10) {
    threshold <- 10
  } else {
    threshold <- round(as.numeric(peaks$cov[1])*.8) 
  }
  
  for (i in seq_len(nrow(plateaus))) {
    matched <- FALSE
    for (j in seq_len(nrow(peaks))) {
      if (abs(plateaus$cov[i] - peaks$cov[j]) <= threshold) {
        matched <- TRUE
        break
      }
    }
    if (!matched) {
      merged <- rbind(merged, plateaus[i, ])
      merged <-merged[order(merged$cov), ]
    }
  }
  
  return(merged)
}

# Function to calculate ploidy
get_ploidy <- function(peaks) {
  # Return an empty character vector if there are no peaks or the first peak's coverage is greater than 120
  if (nrow(peaks) == 0 || as.numeric(peaks$cov[1]) > 120) return(character(0))
  
  n_cov <- peaks$cov[1]  # Get the first coverage value
  
  ploidy <- sapply(peaks$cov, function(x) {
    if (x == n_cov) {
      return("n")  # Assign "n" for the first coverage value
    } else {
      ploidy_level <- round(x / n_cov)  # Calculate ploidy level
      return(paste0(ploidy_level, "n"))  # Append "n" to the calculated ploidy level
    }
  })
  
  return(ploidy)
}

# Search for the file ending with _histo.txt in the current directory
file_pattern <- "_histo.txt$"
data_file <- list.files(pattern = file_pattern)

  
# Assume the total_bases file is in the same directory and has a fixed name
total_bases <- 'total_number_bases.txt$'
total_bases_path<- list.files(pattern = total_bases)
  
species_info <- strsplit(data_file, "_")[[1]]
species_name <-  species_info[1]
kmer <- species_info[2]
  
# Read data
data <- read.table(data_file[1])
total_bases <- read.table(total_bases_path[1])
colnames(data) <- c("cov", "freq")

# Find peaks and plateaus
peaks <- find_peaks(data)
plateaus <- find_plateaus(data)

# Apply the custom merge function
peaks_and_plateaus <- merge_peaks_and_plateaus(peaks, plateaus)
# Generate output file based on number of peaks and plateaus present 
if (nrow(peaks_and_plateaus) == 0) {
    # Write output text file
  writeLines(paste(species_name, "K-mer: ", kmer,"\nNo peaks or plateaus detected"), paste0(species_name,"_K_mers.txt"))
    
  # Create a basic plot without peaks and plateaus pdf
  y_limt <- sort(data$freq, decreasing = TRUE)
  pdf(filename = paste0(species_name, "_K_mers.pdf"), width = 24, height = 14)
  plot(data$cov, data$freq, type = "l", xlim = c(0, 250), ylim = c(0, y_limt[2]),
       xlab = "Coverage", ylab = "Frequency", main = paste(species_name,  "K-mer: ", kmer, "\nNo peaks or plateaus detected"))
  dev.off()
  # Create a basic plot without peaks and plateaus png
  png(filename = paste0(species_name,"_K_mers.png"), width = 24, height = 14, units = "cm", res = 300)  
  plot(data$cov, data$freq, type = "l", xlim = c(0, 250), ylim = c(0, y_limt[2]),
       xlab = "Coverage", ylab = "Frequency", main = paste(species_name,  "K-mer: ", kmer, "\nNo peaks or plateaus detected"))
  dev.off()
} else {
  # Get ploidy levels
  ploidy_levels <- get_ploidy(peaks_and_plateaus)
  if (length(ploidy_levels) == 0) {
    peaks_and_plateaus$ploidy <- NA
    peaks_and_plateaus$label <- NA
  } else {
    peaks_and_plateaus$ploidy <- ploidy_levels
    peaks_and_plateaus$label <- paste(peaks_and_plateaus$ploidy, peaks_and_plateaus$cov, sep = ":  ")
  }
  # Estimate n/2n info
  if (nrow(peaks_and_plateaus) > 1) {
    if (peaks_and_plateaus$ploidy[2]=="2n") {
      length <- as.numeric(total_bases[2, 1]) / as.numeric(peaks_and_plateaus$cov[2])
      cov<-as.numeric(peaks_and_plateaus$cov[2])
      Mb_value <- as.numeric(length) / 1e6
      formatted_value <- format(round(Mb_value, 2), big.mark = " ", decimal.mark = ".", nsmall = 2)
      length_value <- paste(formatted_value, "Mb")
      length_lable <- "Genome Length:"
      assembly_cov <- tail(peaks_and_plateaus, n = 1)
      assembly_length <- as.numeric(total_bases[2, 1]) / as.numeric(assembly_cov$cov)
      assembly_Mb_value <- as.numeric(assembly_length) / 1e6
      assembly_formatted_value <- format(round(assembly_Mb_value, 2), big.mark = " ", decimal.mark = ".", nsmall = 2)
      assembly_length_value <- paste(assembly_formatted_value, "Mb")
    } else {
      length <- as.numeric(total_bases[2, 1]) / round(as.numeric(peaks_and_plateaus$cov[1])*2)
      cov<-round(as.numeric(peaks_and_plateaus$cov[1])*2)
      Mb_value <- as.numeric(length) / 1e6
      formatted_value <- format(round(Mb_value, 2), big.mark = " ", decimal.mark = ".", nsmall = 2)
      length_value <- paste(formatted_value, "Mb")
      length_lable <- "Genome Length:"
      assembly_cov <- tail(peaks_and_plateaus, n = 1)
      assembly_length <- as.numeric(total_bases[2, 1]) / as.numeric(assembly_cov$cov)
      assembly_Mb_value <- as.numeric(assembly_length) / 1e6
      assembly_formatted_value <- format(round(assembly_Mb_value, 2), big.mark = " ", decimal.mark = ".", nsmall = 2)
      assembly_length_value <- paste(assembly_formatted_value, "Mb")
    }
  } else {
    length <- as.numeric(total_bases[2, 1]) /round(as.numeric(peaks_and_plateaus$cov[1])*2)
    cov<-round(as.numeric(peaks_and_plateaus$cov[1])*2)
    Mb_value <- as.numeric(length) / 1e6
    formatted_value <- format(round(Mb_value, 2), big.mark = " ", decimal.mark = ".", nsmall = 2)
    length_value <- paste(formatted_value, "Mb")
    length_lable <- "Haploid Length:"
    assembly_length <- (as.numeric(total_bases[2, 1]) / as.numeric(peaks_and_plateaus$cov[1]))/2
    assembly_Mb_value <- as.numeric(assembly_length) / 1e6
    assembly_formatted_value <- format(round(assembly_Mb_value, 2), big.mark = " ", decimal.mark = ".", nsmall = 2)
    assembly_length_value <- paste(assembly_formatted_value, "Mb")
  }
  # Sort peaks_and_plateaus by frequency in decreasing order
  y_limt <- sort(peaks_and_plateaus$freq, decreasing = TRUE)
  
  # Adjust y_limt
  if (length(y_limt) == 0 || y_limt[1] < 1e6) {
    y_limt <- round(as.numeric(data$freq[4]))
  } else {
    y_limt <- round(as.numeric(y_limt[1]) * 1.1)
  }
  
  # Set x_limt based on the last peak's coverage
  last_peak <- sort(peaks_and_plateaus$cov, decreasing = TRUE)
  if (length(last_peak) > 0 && 250 < as.numeric(last_peak[1]) && as.numeric(last_peak[1]) < 1000) {
    x_limt <- as.numeric(last_peak[1]) + 100
  } else {
    x_limt <- 250
  }
  # Write Output file
  writeLines(paste(species_name,"K-mer:",kmer,"\nEstimated",length_lable,length_value,"\nEstimated Coverage:",cov,"\nExpected Assembly Length:",assembly_length_value), paste0(species_name,"_K_mers.txt"))
  # Create a basic plot with peaks and or plateaus pdf
  pdf(file = paste0(species_name, "_K_mers.pdf"), width = 24, height = 14)
  plot(data$cov, data$freq, type = "l", xlim = c(0, x_limt), ylim = c(0, y_limt),
       xlab = "Coverage", ylab = "Frequency", 
       main = paste(species_name,"K-mer:",kmer,"\nEstimated",length_lable,length_value,"\nExpected Assembly Length:",assembly_length_value))  
  abline(v = peaks_and_plateaus$cov, lty = 2)
  text(x = peaks_and_plateaus$cov, y = peaks_and_plateaus$freq, labels = peaks_and_plateaus$label,
       pos = 3, cex = 0.8)
  dev.off()
  # Create a basic plot with peaks and or plateaus png
  png(filename = paste0(species_name,"_K_mers.png"), width = 24, height = 14, units = "cm", res = 300)
  plot(data$cov, data$freq, type = "l", xlim = c(0, x_limt), ylim = c(0, y_limt),
       xlab = "Coverage", ylab = "Frequency", 
       main = paste(species_name,"K-mer:",kmer,"\nEstimated",length_lable,length_value,"\nExpected Assembly Length:",assembly_length_value))
  abline(v = peaks_and_plateaus$cov, lty = 2)
  text(x = peaks_and_plateaus$cov, y = peaks_and_plateaus$freq, labels = peaks_and_plateaus$label,
       pos = 3, cex = 0.8)
  dev.off()
}
