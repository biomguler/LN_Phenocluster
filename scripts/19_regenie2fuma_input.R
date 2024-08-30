#############################################################################
#Title:    Regenie sumstats to fuma input
#Function: Create input for fuma
#Author:   Murat Guler
#Time:     March 26th 2024
#############################################################################
rm(list=ls())
gc()
#############################################################################
#required packages
library(data.table)
library(tidyverse)


# Working dir
setwd("/your/sumstats/dir")

Files <- list.files(path = "/your/sumstats/dir",
                    pattern = "\\_out_A1FREQ0p001_info0p9.txt$",
                    full.names = TRUE,
                    include.dirs = TRUE,
                    recursive = TRUE)

# List of phenotypes
phenos <- basename(dirname(Files))

# Loop through phenotypes
for (i in 1:length(Files)) {
  # Read the file
  df <- fread(Files[i], sep = " ")
  df <- df %>% select(CHROM, GENPOS, ALLELE0, ALLELE1, BETA, SE, P)
  # Construct output path based on phenotype
  out_path <- file.path(phenos[i], paste0(phenos[i], "_09AF001_fuma.gz"))
  # Save files
  fwrite(df, out_path, sep = " ", compress = "gzip")
}

#End