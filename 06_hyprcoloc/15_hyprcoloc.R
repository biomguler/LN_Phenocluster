#!/usr/bin/env Rscript
#############################################################################
# Title:    HyPrColoc Analysis Using Pre-Merged Data
# Author:   Murat Guler
# Time:     December 3rd 2024
#############################################################################

rm(list = ls())
gc()

#############################################################################
# Required packages
#############################################################################
suppressMessages({
  library(data.table)
  library(tidyverse)
  library(hyprcoloc)
})

#############################################################################
# Parse command-line arguments
#############################################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 02_run_hyprcoloc.R <merged_data_path> <region_path> [pval_threshold] <outname>")
}

# Arguments:
#   merged_data_path: path to the pre-saved merged_data.RData file.
#   region_path: path to region file (CSV/TSV with columns: LOC, CHR, START, STOP)
#   pval_threshold (optional): p-value threshold for filtering traits in the region.
merged_data_path <- args[1]
region_path      <- args[2]
pval_threshold   <- if (length(args) >= 3) as.numeric(args[3]) else 1e-5
outname   <- if (length(args) >= 4) as.character(args[4]) else "out"

message("Using merged data file: ", merged_data_path)
message("Using region file: ", region_path)
message("Using p-value threshold: ", pval_threshold)
message("Using output name: ", outname)
#############################################################################
# Step 1: Load merged data and region file
#############################################################################
load(merged_data_path)  # Loads the object 'merged_data'
if (!exists("merged_data")) {
  stop("merged_data not found in ", merged_data_path)
}
region <- fread(region_path) %>% as.data.frame()

#############################################################################
# Step 2: For each region, subset the data and run HyPrColoc
#############################################################################
results_hyprcoloc <- list()
info_region <- list()

# Check for LSF job array index (if used).
job_index_str <- Sys.getenv("LSB_JOBINDEX")
if (job_index_str != "") {
  region_indices <- as.numeric(job_index_str)
  message("Processing region index: ", region_indices)
} else {
  region_indices <- 1:nrow(region)
}

for (i in region_indices) {
  # Filter the merged data for the current region by CHR and POS.
  merged_region <- merged_data %>% 
    filter(CHR == region$CHR[i],
           POS >= region$START[i],
           POS <= region$STOP[i])
  
  # Determine which traits to keep.
  # For each p-value column (e.g., "trait.PVAL"), keep the trait if at least one SNP
  # in the region has a p-value below the threshold.
  pval_cols <- grep("\\.PVAL$", colnames(merged_region), value = TRUE)
  keep_traits <- c()
  for (pcol in pval_cols) {
    if (any(merged_region[[pcol]] < pval_threshold)) {
      trait_name <- sub("\\.PVAL$", "", pcol)
      keep_traits <- c(keep_traits, trait_name)
    }
  }
  
  # Skip the region if no traits remain.
  if (length(keep_traits) < 1) next
  
  # Select only the columns for the traits that passed the p-value filter.
  beta_columns <- paste(keep_traits, "BETA", sep = ".")
  se_columns   <- paste(keep_traits, "SE", sep = ".")
  
  # Create the BETA matrix.
  beta_matrix <- as.matrix(merged_region[, ..beta_columns])
  rownames(beta_matrix) <- merged_region$SNP
  colnames(beta_matrix) <- sub("\\.BETA$", "", colnames(beta_matrix))
  
  # Create the SE matrix.
  se_matrix <- as.matrix(merged_region[, ..se_columns])
  rownames(se_matrix) <- merged_region$SNP
  colnames(se_matrix) <- sub("\\.SE$", "", colnames(se_matrix))
  # create region name
  region_name <- paste0("Region_", i)
  # (Optional) Collect region-level info.
  info <- data.frame(Region = i,
                     nSNP = nrow(beta_matrix),
                     Traits = paste(colnames(beta_matrix), collapse = ", "))
  info_region[[region_name]] <- info
  
  # Run HyPrColoc.
  betas  <- beta_matrix
  ses    <- se_matrix
  traits <- colnames(beta_matrix)
  rsid   <- rownames(beta_matrix)
  
  res <- hyprcoloc(betas, ses, trait.names = traits, 
                   snp.id = rsid, snpscores = TRUE, bb.alg = TRUE, 
                   bb.selection = "alignment")
  results_hyprcoloc[[region_name]] <- res
}

#############################################################################
# Step 3: Save the results and region info for further use
#############################################################################
if (length(results_hyprcoloc) > 0) {
  if (length(region_indices) == 1) {
    outfile <- paste0(outname, "_", sprintf("hyprcoloc_results_region_%d.RData", region_indices))
  } else {
    outfile <- paste0(outname, "_", "hyprcoloc_test_results.RData")
  }
  save(results_hyprcoloc, info_region, file = outfile)
  message("Results saved to: ", outfile)
} else {
  message("No results produced; nothing to save.")
}