#############################################################################
# Title:    HyPrColoc Results Aggregation Script
# Author:   Murat Güler
# Date:     December 3rd 2024
#############################################################################

rm(list = ls())
gc()

#############################################################################
# Required packages
#############################################################################
library(data.table)

#############################################################################
# Initialize containers for all result parts
#############################################################################
all_PIP <- list()
all_info <- list()
all_results <- list()

# Loop through all region result files
for (i in 1:321) {
  region_file <- paste0("LN_regions_hyprcoloc_results_region_", i, ".RData")
  
  if (!file.exists(region_file)) {
    warning("File not found: ", region_file)
    next
  }
  
  load(region_file)  # Loads: results_hyprcoloc and info_region

  # Extract region info
  all_info[[i]] <- as.data.frame(info_region[[1]])
  
  # Check if PIP results exist
  snpscores <- results_hyprcoloc[[1]][["snpscores"]]
  if (length(snpscores) == 0) next
  
  # Flatten PIP list and add metadata
  pip_df <- data.frame(
    SNP = names(unlist(snpscores)),
    PIP = unlist(snpscores),
    Region = i
  )
  all_PIP[[i]] <- pip_df
  
  # Extract coloc results
  coloc_res <- as.data.frame(results_hyprcoloc[[1]][["results"]])
  coloc_res$Region <- i
  all_results[[i]()]()_
