#!/usr/bin/env Rscript
#############################################################################
# Title:    flames_index files generator
# Author:   Murat Guler
# Date:     December 3rd, 2024
#############################################################################
# Notes:
# 1: If there is not cs r2 > 0.8 && Pgwas < 1e-5 taken

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

# Read the file
data <- fread("indexfile.txt")

# Get unique values of locus_set
unique_pheno <- unique(data$pheno)

# Loop through each unique value and write separate files
for (set in unique_pheno) {
  subset_data <- filter(data, pheno == set) 
  subset_data <- subset_data %>% 
    mutate(GenomicLocus = c(1:nrow(subset_data))) %>% 
    select(Filename,	GenomicLocus,	Annotfiles)
  out_file <- paste0(set, "_indexfile.txt")
  fwrite(subset_data,  out_file, sep = "\t")
}

print("Files have been saved separately for each locus_set.")
