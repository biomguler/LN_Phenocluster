#!/usr/bin/env Rscript
#############################################################################
# Title:    susie output merger
# Author:   Murat Guler
# Date:     December 3rd, 2024
#############################################################################

rm(list = ls())  # Clear environment
gc()  # Run garbage collection

#############################################################################
# Required packages
#############################################################################
suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
})

#############################################################################
# Read the list of finemap outputs
#############################################################################
finemap_out_list <- fread("finemap_out.list")

# Initialize lists to store data
credible_sets <- list()
credible_info <- list()

#############################################################################
# Process Credible Sets
#############################################################################
for (i in 1:nrow(finemap_out_list)) {
  file <- as.character(finemap_out_list[i, 3][[1]])  # Extract filename
  
  if (!file.exists(file)) {
    warning(paste("Skipping", file, "- file does not exist"))
    next
  }
  
  df <- fread(file)
  
  # Ensure 'cs' column exists before filtering
  if (!"cs" %in% colnames(df)) {
    warning(paste("Skipping", file, "- column 'cs' not found"))
    next
  }
  
  # Filter out rows where cs == -1
  filtered_df <- df %>%
    filter(cs != -1) %>%
    mutate(locus = as.character(finemap_out_list[i, 4][[1]]))
  
  # If filtering removes all rows, use the original (unfiltered) df
  if (nrow(filtered_df) > 0) {
    credible_sets[[i]] <- filtered_df
  } else {
    credible_sets[[i]] <- df  %>%
      mutate(locus = as.character(finemap_out_list[i, 4][[1]]))
    warning(paste("Using unfiltered data for", file, "- no rows remained after filtering"))
  }
}

# Combine all credible sets into one data frame
final_credible_sets <- rbindlist(credible_sets, use.names = TRUE, fill = TRUE, idcol = "source")


#############################################################################
# Process Credible Info
#############################################################################
for (i in 1:nrow(finemap_out_list)) {
  file <- as.character(finemap_out_list[i, 5][[1]])  # Extract filename
  
  if (!file.exists(file)) {
    warning(paste("Skipping", file, "- file does not exist"))
    next
  }
  
  df <- fread(file)
  
  # Ensure it has data and add locus column
  df <- df %>%
    mutate(locus = as.character(finemap_out_list[i, 4][[1]]))
  
  credible_info[[i]] <- df
}

# Combine all credible info into one data frame
final_credible_info <- rbindlist(credible_info, use.names = TRUE, fill = TRUE, idcol = "source")


#############################################################################
# Save Results
#############################################################################
fwrite(final_credible_sets, "final_credible_sets.tab", sep = "\t")
fwrite(final_credible_info, "final_credible_info.tab", sep = "\t")

cat("✅ Successfully saved final_credible_sets.tab and final_credible_info.tab\n")
