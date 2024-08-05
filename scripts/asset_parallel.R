#############################################################################
# Title:    ASSET parallel
# Function: ASSET analysis
# Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
# Date:     Jan 23th 2022
# Note:     Please let me know if you have any trouble
#############################################################################
# Free R memory and remove prior environment
rm(list=ls())
gc()
#############################################################################
# Packages

# Install required packages, if not installed

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# List of packages to install and load
packages <- c("data.table", "ASSET", "dplyr", "purr", "parallel")

# Install and load each package
for (package in packages) {
  install_if_missing(package)
}

# Load required packages
library(ASSET)
library(data.table)
library(dplyr)
library(parallel)
library(purrr)

# setwd
setwd("/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/")

# Define phenotypes
phenotypes <- c("CLL", "DLBCL", "FL", "HL", "LPL_WM", "MGUS", "MM", "MZL")

# Get data dir
main_dir <- "/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS"

# Initialize an empty list to store paths
paths_list <- list()

# Loop through phenotypes and store paths
for (pheno in phenotypes) {
  path <- paste0(main_dir, "/", pheno)
  paths_list[[pheno]] <- path
}


# Read GWAS sumstat, extract BETA and SE

gwas_list <- list()
for (i in 1:8) {
  pheno <- phenotypes[i]
  sumstats_path <- paste0(paths_list[i], "/", phenotypes[i], "_out_A1FREQ0p001_info0p9.txt")
  sumstats <- fread(sumstats_path, header = TRUE, na.strings="NA", select= c(19,14,15), sep = " ")
  names(sumstats) <- c("SNP", paste0(phenotypes[i], ".BETA"), paste0(phenotypes[i], ".SE"))
  gwas_list[[pheno]] <- sumstats}
 


# Merge sumstats

gwas_merge <- reduce(gwas_list, inner_join, by = "SNP")

# Define the input arguments to h.traits
traits.lab <- c("CLL","DLBCL", "FL","HL","LPL_WM", "MGUS", "MM", "MZL")
N11 <- paste0(main_dir, "/", "asset", "/", "N11.txt")
N10 <- paste0(main_dir, "/", "asset", "/", "N10.txt")
N00 <- paste0(main_dir, "/", "asset", "/", "N00.txt")

N11 <- as.matrix(fread(N11, drop = 1))
N10 <- as.matrix(fread(N10, drop = 1))
N00 <- as.matrix(fread(N00, drop = 1))
rownames(N11) <- traits.lab
rownames(N10) <- traits.lab
rownames(N00) <- traits.lab
cor <- list(N11=N11, N00=N00, N10=N10)
ncase <- diag(N11)
ncntl <- diag(N00)

snps <- as.vector(gwas_merge[, "SNP"])

# Extract column names that end with ".BETA"
beta_columns <- grep("\\.BETA$", names(gwas_merge), value = TRUE)

# Extract the corresponding columns and convert to matrix
beta.hat <- as.matrix(gwas_merge[, beta_columns, with = FALSE])

rownames(beta.hat) <- seq_len(nrow(gwas_merge))

# Extract column names that end with ".BETA"
sigma_columns <- grep("\\.SE$", names(gwas_merge), value = TRUE)

# Extract the corresponding columns and convert to matrix
sigma.hat <- as.matrix(gwas_merge[, sigma_columns, with = FALSE])
rownames(sigma.hat) <- seq_len(nrow(gwas_merge))

# Number of cores to use
num_cores <- 64

# Initialize a cluster
cl <- makeCluster(num_cores)

# Parallel loop using mclapply
summary_list <- mclapply(1:nrow(snps), function(i) {
  res <- h.traits(snps[i], traits.lab, beta.hat[i, , drop = FALSE], sigma.hat[i, , drop = FALSE], ncase, ncntl, cor = cor, cor.numr = FALSE, side = 2, meta = TRUE, zmax.args = NULL, meth.pval = "DLM")
  
  # Get individual data frames
  meta_df <- h.summary(res)[["Meta"]]
  colnames(meta_df)[2:5] <- paste0(colnames(meta_df)[2:5], "_meta")
  subset_1sided_df <- h.summary(res)[["Subset.1sided"]]
  colnames(subset_1sided_df)[2:6] <- paste0(colnames(subset_1sided_df)[2:6], "_1sided")
  subset_2sided_df <- h.summary(res)[["Subset.2sided"]]
  colnames(subset_2sided_df)[2:12] <- paste0(colnames(subset_2sided_df)[2:12], "_2sided")
  
  # Merge the three data frames by the "SNP" column
  merged_df <- merge(meta_df, subset_1sided_df, by = "SNP")
  merged_df <- merge(merged_df, subset_2sided_df, by = "SNP")
  
  # Return the merged data frame
  return(merged_df)
}, mc.cores = num_cores)

# Stop the cluster
stopCluster(cl)

# Combine the data frames within the summary_list into a single data frame
combined_df <- bind_rows(summary_list)

#Save
setwd("/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/asset/")

# Write combined_df to a tab-delimited text file without row names
fwrite(combined_df, "asset_sumstat.txt", sep = "\t", quote = FALSE, row.names = FALSE)



