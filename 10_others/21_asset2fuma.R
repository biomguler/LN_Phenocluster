#############################################################################
#Title:    ASSET sumstats to fuma input
#Function: Create input for fuma from ASSET results
#Author:   Murat Guler
#Time:     March 26th 2024
#############################################################################
rm(list=ls())
gc()
#############################################################################
#required packages
library(data.table)
library(tidyverse)
setwd("/your/asset/sumstats/dir")

# Load ASSET results
ukbb_asset <- fread("asset_sumstat.txt", sep = "\t")


# Filter single phenotype based assocations

ukbb_asset_1sided_raw <- ukbb_asset %>% filter(str_detect(Pheno_1sided, ",")) %>%
  separate(col = SNP, into = c("CHROM", "GENPOS", "REF", "ALT"), sep = ":", remove = FALSE)

ukbb_asset_2sided_raw <- ukbb_asset %>% filter(str_detect(Pheno.1_2sided, ",") | str_detect(Pheno.2_2sided, ",")) %>%
  separate(col = SNP, into = c("CHROM", "GENPOS", "REF", "ALT"), sep = ":", remove = FALSE)

# create fuma summary files

ukbb_asset_1sided_raw <- ukbb_asset_1sided_raw %>% 
  mutate (BETA = log(as.numeric(OR_1sided)), SE = (log(as.numeric(CI.high_1sided)) - log(as.numeric(CI.low_1sided))) / 3.92) %>% 
  select(CHROM, GENPOS, REF, ALT, BETA, SE, Pvalue_1sided) 

ukbb_asset_1sided_raw <- ukbb_asset_1sided_raw %>%
  filter(is.finite(BETA) & is.finite(SE))

fwrite(ukbb_asset_1sided_raw, "raw_ukbb_asset_1sided_fuma.gz", sep = " ", compress = "gzip")

ukbb_asset_2sided_raw <- ukbb_asset_2sided_raw %>% 
  select(CHROM, GENPOS, REF, ALT, Pvalue_2sided)

fwrite(ukbb_asset_2sided_raw, "raw_ukbb_asset_2sided_fuma.gz", sep = " ", compress = "gzip")
