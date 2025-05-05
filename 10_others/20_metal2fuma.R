#############################################################################
#Title:    Metal sumstats to fuma input
#Function: Create input for fuma from metal results
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
setwd("/your/metalsumstats/dir")

input_paths <- list.files(path = "/your/metalsumstats/dir",
                          pattern = "\\.tbl$",
                          full.names = TRUE,
                          include.dirs = TRUE,
                          recursive = TRUE)

output_names <- list.files(pattern="\\.tbl$")
output_names <- gsub("\\.tbl$", "", output_names)



# Loop through phenotypes
for (i in 1:length(input_paths)) {
  # Read the file
  df <- fread(input_paths[i], sep = "\t")
  df <- df %>% filter(HetDf >= mean(nchar(df$Direction))-1, HetPVal > 0.05, HetISq < 80) %>%
    arrange(`P-value`)
  # chr and position, select
  df <- df %>%
    separate(col = MarkerName, into = c("CHROM", "GENPOS", "REF", "ALT"), sep = ":", remove = FALSE) %>%
    select(CHROM, GENPOS, Allele1, Allele2, Effect, StdErr, `P-value`) %>%
    rename(EA = Allele1, NEA= Allele2, Beta=Effect, SE = StdErr, P = `P-value`)
  # Construct output name
  out_path <- paste0(output_names[i], "_fuma.gz")
  # Save files
  fwrite(df, out_path, sep = " ", compress = "gzip")
}



