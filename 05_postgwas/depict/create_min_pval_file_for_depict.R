# makes 2 text files for depict 
library(dplyr)
library(tidyr)
library(data.table)
library(openxlsx)
library(readr)
library(stringr)
library(gdata)

setwd("/mnt/polyomica/projects/GPD/discovery_meta/")
path = "/mnt/polyomica/projects/GPD/PostGWAS_v2/"

thr1 = 1e-5/12
thr2 = 5e-8/12

# Define the IGP numbers you're using
igp_numbers <- c(2,3,4,5,7,8,9,10,11,12,15,20)
list_files <- paste0("IGP", igp_numbers, "/igp", igp_numbers, ".csv")
names <- c('chr', 'bp', 'rs_id', 'p_ma')

# Initialize empty list to store data
data_list <- list()

# First verify all files exist
missing_files <- list_files[!file.exists(list_files)]
if (length(missing_files) > 0) {
  warning(paste("The following files are missing:", paste(missing_files, collapse=", ")))
}

# Process each file with error handling
for (i in seq_along(igp_numbers)) {
  current_file <- list_files[i]
  current_igp <- igp_numbers[i]
  
  cat("Processing", current_file, "...\n")
  
  tryCatch({
    # Read the file
    temp_data <- fread(current_file, select = names)
    
    # Filter by p-value threshold
    temp_data <- temp_data[temp_data$p_ma < thr1, ]
    
    if (nrow(temp_data) > 0) {
      # Rename p-value column
      names(temp_data)[4] <- paste0("p_ma", current_igp)
      data_list[[as.character(current_igp)]] <- temp_data
    } else {
      cat("No significant SNPs in", current_file, "\n")
    }
  }, error = function(e) {
    warning(paste("Failed to process", current_file, ":", e$message))
  })
}

# Combine all data frames
if (length(data_list) > 0) {
  # Start with first dataset
  data <- data_list[[1]]
  
  # Merge remaining datasets
  if (length(data_list) > 1) {
    for (i in 2:length(data_list)) {
      data <- merge(data, data_list[[i]], by = names[1:3], all = TRUE)
    }
  }
  
  # Calculate minimum p-value
  data <- as.data.frame(data)
  data$p <- apply(data[, startsWith(names(data), "p_")], 1, min, na.rm = TRUE)
  
  # Prepare output
  data_small <- data[, c(3, 1, 2, ncol(data))]
  data_small$p <- signif(data_small$p, 3)
  data_very_small <- data_small[data_small$p < thr2, ]
  
  # Write output files
  write.table(data_small, paste0(path, "for_depict_all_e-5.txt"), 
              sep = ",", quote = FALSE, row.names = FALSE) 
  write.table(data_very_small, paste0(path, "for_depict_all_5e-8.txt"), 
              sep = ",", quote = FALSE, row.names = FALSE)
  
  cat("Successfully created output files.\n")
} else {
  warning("No data was processed - check input files and thresholds")
}
