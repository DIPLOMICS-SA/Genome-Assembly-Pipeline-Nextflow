##################################################
#R script for K-mer analysis
#Wilku Meyer 
#contact: wilku@cengen.ac.za
##################################################

# Function to find peaks that are more than 80% units of the first peak apart
# and consider 5 values for comparison
find_peaks <- function(data) {
  peaks <- data.frame(cov = numeric(), freq = numeric())
  for (i in 4:(nrow(data) - 3)) {
    if (data[i, "cov"] < 1000 &&
        data[i, "freq"] > 1500000 &&
        data[i, "freq"] > data[i - 2, "freq"] &&
        data[i, "freq"] > data[i - 1, "freq"] &&
        data[i, "freq"] > data[i + 1, "freq"] &&
        data[i, "freq"] > data[i + 2, "freq"] &&
        data[i, "freq"] > data[i + 3, "freq"]
    ) {
      if (nrow(peaks) == 0 || data[i, "cov"] > 0.8 * tail(peaks$cov, 1)) {
        peaks <- rbind(peaks, data[i, ])
      }
    }
  }
  return(peaks)
}

# Function to find plateaus that are more than 10 units apart
# and are 4 values or longer
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
    while (end < n && data[end, "freq"] > 2500000) {
      diff <- (data[end, "freq"] - data[end + 1, "freq"]) / data[end, "freq"]
      if (diff > -0.01 && diff < 0.01) {
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
      plateaus <- rbind(plateaus, data[round((start + end) / 2), , drop = FALSE])
      plateaus_end <- rbind(plateaus_end, data[end, , drop = FALSE])
    }
    
    # Move to the next potential plateau
    i <- end + 1
  }
  
  return(plateaus)
}

# Custom merge function
merge_peaks_and_plateaus <- function(peaks, plateaus) {
  if (nrow(peaks) == 0) {
    return(plateaus)
  }
  if (nrow(plateaus) == 0) {
    return(peaks)
  }

  # Determine threshold
  if (is.na(peaks$cov[1]) || peaks$cov[1] < 10) {
    threshold <- 10
  } else {
    threshold <- round(as.numeric(peaks$cov[1]) * 0.8)
  }

  # Filter plateaus that are not within threshold of any peak
  unmatched_plateaus <- plateaus[
    sapply(plateaus$cov, function(pc) {
      all(abs(pc - peaks$cov) > threshold)
    }),
    ,
    drop = FALSE
  ]

  # Combine and sort
  merged <- rbind(peaks, unmatched_plateaus)
  merged <- merged[order(merged$cov), , drop = FALSE]

  return(merged)
}

# Function to calculate ploidy
get_ploidy <- function(peaks) {
  # Check if peaks is empty or if the first coverage value is greater than 120
  if (nrow(peaks) == 0 || as.numeric(peaks$cov[1]) > 120) {
    return(data.frame(cov = numeric(0), ploidy = numeric(0)))
  }
  
  # Handle the case where there is only one peak
  if (nrow(peaks) == 1) {
    return(data.frame(cov = peaks$cov[1], freq = peaks$freq[1], ploidy = "?n",
                      label = paste0("?n : ", peaks$cov[1])))
  }
  
  complete_cov <- peaks$cov[1]  # First coverage value
  peak_multiple <- round(peaks$cov / complete_cov)
  
  lowest_cov <- peaks$cov[1]
  lowest_freq <- peaks$freq[1]
  lowest_multiple <- peak_multiple[1]
  
  # Second coverage
  second_cov <- peaks$cov[2]
  second_freq <- peaks$freq[2]
  second_multiple <- peak_multiple[2]
  
  # If only two peaks
  if (nrow(peaks) == 2) {
    return(data.frame(
      cov = c(lowest_cov, second_cov),
      freq = c(lowest_freq, second_freq),
      ploidy = c(paste0(lowest_multiple, "n"), paste0(second_multiple, "n")),
      label = c(paste0(lowest_multiple, "n : ", lowest_cov),
                paste0(second_multiple, "n : ", second_cov))
    ))
  }
  
  # If more than two peaks, include the last
  highest_cov <- peaks$cov[nrow(peaks)]
  highest_freq <- peaks$freq[nrow(peaks)]
  highest_multiple <- peak_multiple[nrow(peaks)]
  
  return(data.frame(
    cov = c(lowest_cov, second_cov, highest_cov),
    freq = c(lowest_freq, second_freq, highest_freq),
    ploidy = c("1n", paste0(second_multiple, "n"), paste0(highest_multiple, "n")),
    label = c(
      paste0("1n : ", lowest_cov),
      paste0(second_multiple, "n : ", second_cov),
      paste0(highest_multiple, "n : ", highest_cov)
    )
  ))
}

# ----------------- MAIN SCRIPT ------------------

# Search for the file ending with _histo.txt in the current directory
file_pattern <- "_histo.txt$"
data_file <- list.files(pattern = file_pattern)

# total_bases file is in the same directory and has a fixed name from KMC script
total_bases <- readLines("total_number_bases.txt")
species_info <- strsplit(data_file, "_")[[1]]
species_name <- paste(species_info[1:2], collapse = "_")
kmer <- species_info[3]

# Read data
data <- read.table(data_file[1])
colnames(data) <- c("cov", "freq")
# Keep only freq > 0
data <- data[data$freq > 0, ]

# Find peaks & plateaus, then merge
peaks <- find_peaks(data)
plateaus <- find_plateaus(data)
peaks_and_plateaus <- merge_peaks_and_plateaus(peaks, plateaus)

# --------------- LOGIC: NO PEAKS vs PEAKS ---------------

if (nrow(peaks_and_plateaus) == 0) {
  # Case: No peaks or plateaus
  writeLines(paste0(
    species_name, ": k-mer= ", kmer,
    "\nTotal input bases ", total_bases[2],
    "\nPeaks or Plateaus detected= 0",
    "\nPloidy= ?",
    "\nSize= ? at ? Coverage\n"
  ), paste0("k_mers_Stats_", species_name, ".txt"))
  
  # Basic PDF
  y_limt <- sort(data$freq, decreasing = TRUE)
  pdf(file = paste0(species_name, "_k_mers.pdf"), width = 24, height = 14)
  plot(data$cov, data$freq, type = "l",
       xlim = c(0, 250),
       ylim = c(0, y_limt[2]),
       xlab = "Coverage", ylab = "Frequency",
       main = paste(species_name,
                    "\nK-mer: ", kmer,
                    "Total input bases ", total_bases[2],
                    "\nNo peaks or plateaus detected"))
  dev.off()
  
  # Basic PNG
  png(filename = paste0(species_name, "_k_mers.png"),
      width = 24, height = 14, units = "cm", res = 300)
  plot(data$cov, data$freq, type = "l",
       xlim = c(0, 250),
       ylim = c(0, y_limt[2]),
       xlab = "Coverage", ylab = "Frequency",
       main = paste(species_name,
                    "\nk-mer: ", kmer,
                    "Total input bases ", total_bases[2],
                    "\nNo peaks or plateaus detected"))
  dev.off()
  
  # -----------------NEW ASCII -----------------
  file_conn <- file(paste0("k_mers_Stats_", species_name, ".txt"), open = "a")
  writeLines("\nASCII Plot Representation:\n", file_conn)
  
  # 1) Choose coverage window (just a simple default range if no peaks):
  ascii_start <- 1
  ascii_end <- 100
  
  top_rows <- data[data$cov >= ascii_start & data$cov <= ascii_end, ]
  
  
  # 2) Filter out top outliers in freq
  Q3 <- quantile(top_rows$freq, 0.75)
  IQR_val <- IQR(top_rows$freq)
  upper_bound <- Q3 + 1.5 * IQR_val
  filtered_freqs <- top_rows$freq[top_rows$freq <= upper_bound]
  ave_freq <- if (length(filtered_freqs) > 0) mean(filtered_freqs) else 1
  
  # 3) Print lines with scaled stars
  for (i in seq_len(nrow(top_rows))) {
    freq_val <- top_rows$freq[i]
    cov_val  <- top_rows$cov[i]
    scaled_count <- round((freq_val / ave_freq) * 60)
    stars <- paste0(rep("*", min(scaled_count, 60)), collapse = "")
    marker <- if (scaled_count > 60) ">" else ""
    line <- paste0(cov_val, stars, marker, freq_val)
    writeLines(line, file_conn)
  }
  
  close(file_conn)
  
} else {
  # Case: We have peaks or plateaus
  ploidy_levels <- get_ploidy(peaks_and_plateaus)
  
  # 1) Estimate n/2n, etc.
  if (nrow(peaks_and_plateaus) > 1) {
    if (ploidy_levels$ploidy[2] == "2n") {
      length <- as.numeric(total_bases[2]) / as.numeric(peaks_and_plateaus$cov[2])
      cov    <- as.numeric(peaks_and_plateaus$cov[2])
      ploidy_label <- "Ploidy= 2n ="
    } else {
      length <- as.numeric(total_bases[2]) / round(as.numeric(peaks_and_plateaus$cov[1]) * 2)
      cov    <- round(as.numeric(peaks_and_plateaus$cov[1]) * 2)
      ploidy_label <- "Ploidy= 2n ="
    }
  } else {
    length <- as.numeric(total_bases[2]) / round(as.numeric(peaks_and_plateaus$cov[1]))
    cov    <- as.numeric(peaks_and_plateaus$cov[1])
    ploidy_label <- "Ploidy="
  }
  
  Gb_value <- length / 1e9
  formatted_value <- format(round(Gb_value, 2), big.mark = " ",
                            decimal.mark = ".", nsmall = 2)
  length_value <- paste0(formatted_value, " Gb")
  
  assembly_cov <- tail(peaks_and_plateaus, 1)
  assembly_length <- as.numeric(total_bases[2]) / as.numeric(assembly_cov$cov)
  assembly_Gb_value <- assembly_length / 1e9
  assembly_formatted_value <- format(round(assembly_Gb_value, 2),
                                     big.mark = " ", decimal.mark = ".", nsmall = 2)
  assembly_length_value <- paste0(assembly_formatted_value, " Gb")
  
  ploidy <- if (nrow(ploidy_levels) >= 2) ploidy_levels$ploidy[2] else ploidy_levels$ploidy[1]
  
  # 2) Sort peaks_and_plateaus by freq descending for y-limt
  y_limt <- sort(peaks_and_plateaus$freq, decreasing = TRUE)
  if (length(y_limt) == 0 || y_limt[1] < 1e6) {
    y_limt <- round(as.numeric(data$freq[4]))
  } else {
    y_limt <- round(as.numeric(y_limt[1]) * 1.1)
  }
  
  # 3) Coverage limit for PDF/PNG
  last_peak <- sort(peaks_and_plateaus$cov, decreasing = TRUE)
  if (length(last_peak) > 0 && 250 < as.numeric(last_peak[1]) &&
      as.numeric(last_peak[1]) < 1000) {
    x_limt <- as.numeric(last_peak[1]) + 100
  } else {
    x_limt <- 250
  }
  
  # 4) Write final stats to text
  writeLines(
    paste0(
      species_name, ": k-mer=", kmer,
      "\nTotal input bases ", total_bases[2],
      "\nPeaks or Plateaus detected=", nrow(peaks_and_plateaus),
      "\n", ploidy_label, ploidy,
      "\n2n Genome Length=", length_value, " at ", cov, " X Coverage",
      "\nExpected Assembly Length if fully collapsed=", assembly_length_value,
      " at ", assembly_cov$cov, " X Coverage"
    ),
    paste0("k_mers_Stats_", species_name, ".txt")
  )
  
  # 5) Now generate ASCII with the **new** approach (replacing original).
  file_conn <- file(paste0("k_mers_Stats_", species_name, ".txt"), open = "a")
  writeLines("\nASCII Plot Representation:\n", file_conn)
  
  # coverage window
  ascii_start <- round(head(peaks_and_plateaus$cov, 1) * 0.7)
  ascii_end   <- round(tail(peaks_and_plateaus$cov, 1) * 1.1)
  
  if (ascii_end < 100) {
    top_rows <- head(data, 100)
  } else {
    top_rows <- data[data$cov >= ascii_start & data$cov <= ascii_end, ]
  }
  
  # Filter out top freq outliers
  if (nrow(top_rows) > 0) {
    Q3 <- quantile(top_rows$freq, 0.75)
    IQR_val <- IQR(top_rows$freq)
    upper_bound <- Q3 + 1.5 * IQR_val
    filtered_freqs <- top_rows$freq[top_rows$freq <= upper_bound]
    ave_freq <- if (length(filtered_freqs) > 0) mean(filtered_freqs) else 1
    
    for (i in seq_len(nrow(top_rows))) {
      freq_val <- top_rows$freq[i]
      cov_val  <- top_rows$cov[i]
      scaled_count <- round((freq_val / ave_freq) * 60)
      stars <- paste(rep("*", min(scaled_count, 60)), collapse = "")
      marker <- if (scaled_count > 60) ">" else ""
      line <- paste0(cov_val, stars, marker, freq_val)
      writeLines(line, file_conn)
    }
  } else {
    writeLines("(No data in chosen coverage window.)", file_conn)
  }
  close(file_conn)
  
  # 6) Create a basic plot with peaks/plateaus (PDF)
  pdf(file = paste0(species_name, "_k_mers.pdf"), width = 24, height = 14)
  plot(data$cov, data$freq, type = "l", xlim = c(0, x_limt), ylim = c(0, y_limt),
       xlab = "Coverage", ylab = "Frequency",
       main = paste(
         species_name, ": k-mer=", kmer,
         "Total input bases=", total_bases[2],
         "\nEstimated 2n Genome Length=", length_value,
         "\nExpected Assembly Length=", assembly_length_value
       )
  )
  abline(v = peaks_and_plateaus$cov, lty = 2)
  if (nrow(ploidy_levels) > 0) {
    text(
      x = ploidy_levels$cov,
      y = ploidy_levels$freq,
      labels = ploidy_levels$label,
      pos = 3,
      cex = 0.8
    )
  }
  dev.off()
  
  # 7) Create PNG with peaks/plateaus
  png(filename = paste0(species_name, "_k_mers.png"),
      width = 24, height = 14, units = "cm", res = 300)
  plot(data$cov, data$freq, type = "l", xlim = c(0, x_limt), ylim = c(0, y_limt),
       xlab = "Coverage", ylab = "Frequency",
       main = paste(
         species_name, ": k-mer=", kmer,
         "Total input bases=", total_bases[2],
         "\nEstimated 2n Genome Length=", length_value,
         "\nExpected Assembly Length=", assembly_length_value
       )
  )
  abline(v = peaks_and_plateaus$cov, lty = 2)
  if (nrow(ploidy_levels) > 0) {
    text(
      x = ploidy_levels$cov,
      y = ploidy_levels$freq,
      labels = ploidy_levels$label,
      pos = 3,
      cex = 0.8
    )
  }
  dev.off()
}