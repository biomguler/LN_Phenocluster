#!/usr/bin/env Rscript
#############################################################################
# Title:    susie2flames
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
data <- fread("susie_locus.tab")

# Get unique values of locus_set
unique_sets <- unique(data$locus_set)

# Loop through each unique value and write separate files
for (set in unique_sets) {
  subset_data <- filter(data, locus_set == set) 
  subset_data <- subset_data %>% mutate(index = c(1:nrow(subset_data))) %>%
    select (index,	cred1,	prob1)
  out_file <- paste0(set, ".cred")
  fwrite(subset_data,  out_file, sep = "\t")
}

print("Files have been saved separately for each locus_set.")
