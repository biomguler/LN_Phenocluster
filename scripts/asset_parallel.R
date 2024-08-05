#############################################################################
#Title:    asset
#Function: asset analysis
#Author:   Murat Guler
#Time:     Nov 8th 2023
#############################################################################
rm(list=ls())
gc()
#############################################################################

#required packages
library(ASSET)
library(data.table)
library(dplyr)
library(parallel)
library(purrr)
# setwd
setwd("/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/")
# Get data dir
main_dir <- "/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS"
dir_list <- list.dirs(main_dir,recursive = FALSE)
dir_list <- dir_list[dir_list != "/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/LM_files"]
dir_list <- dir_list[dir_list != "/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/asset"]
dir_list <- dir_list[dir_list != "/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/sumResults"]
dir_names <- list.dirs(main_dir,full.names = FALSE, recursive = FALSE)
dir_names <- dir_names[dir_names != "LM_files"]
dir_names <- dir_names[dir_names != "sumResults"]
dir_names <- dir_names[dir_names != "asset"]
df <- data.frame(path=c(dir_list), name=c(dir_names))


# Read GWAS sumstat, extract BETA and SE
CLL <- paste0(df[3, 1], "/", df[3, 2], "_out_A1FREQ0p001_info0p9.txt")
gwas_CLL <- fread(CLL, header = TRUE, na.strings="NA", sep = " ")
gwas_CLL <- gwas_CLL[,c(19,14,15)]
names(gwas_CLL) <- c("SNP", "CLL.BETA", "CLL.SE")

DLBCL <- paste0(df[4, 1], "/", df[4, 2], "_out_A1FREQ0p001_info0p9.txt")
gwas_DLBCL <- fread(DLBCL, header = TRUE, na.strings="NA", sep = " ")
gwas_DLBCL <- gwas_DLBCL[,c(19,14,15)]
names(gwas_DLBCL) <- c("SNP", "DLBCL.BETA", "DLBCL.SE")

FL <- paste0(df[6, 1], "/", df[6, 2], "_out_A1FREQ0p001_info0p9.txt")
gwas_FL <- fread(FL, header = TRUE, na.strings="NA", sep = " ")
gwas_FL <- gwas_FL[,c(19,14,15)]
names(gwas_FL) <- c("SNP", "FL.BETA", "FL.SE")

HL <- paste0(df[7, 1], "/", df[7, 2], "_out_A1FREQ0p001_info0p9.txt")
gwas_HL <- fread(HL, header = TRUE, na.strings="NA", sep = " ")
gwas_HL <- gwas_HL[,c(19,14,15)]
names(gwas_HL) <- c("SNP", "HL.BETA", "HL.SE")

LPL_WM <- paste0(df[9, 1], "/", df[9, 2], "_out_A1FREQ0p001_info0p9.txt")
gwas_LPL_WM <- fread(LPL_WM, header = TRUE, na.strings="NA", sep = " ")
gwas_LPL_WM <- gwas_LPL_WM[,c(19,14,15)]
names(gwas_LPL_WM) <- c("SNP", "LPL_WM.BETA", "LPL_WM.SE")

MGUS <- paste0(df[10, 1], "/", df[10, 2], "_out_A1FREQ0p001_info0p9.txt")
gwas_MGUS <- fread(MGUS, header = TRUE, na.strings="NA", sep = " ")
gwas_MGUS <- gwas_MGUS[,c(19,14,15)]
names(gwas_MGUS) <- c("SNP", "MGUS.BETA", "MGUS.SE")

MM <- paste0(df[11, 1], "/", df[11, 2], "_out_A1FREQ0p001_info0p9.txt")
gwas_MM <- fread(MM, header = TRUE, na.strings="NA", sep = " ")
gwas_MM <- gwas_MM[,c(19,14,15)]
names(gwas_MM) <- c("SNP", "MM.BETA", "MM.SE")

MZL <- paste0(df[13, 1], "/", df[13, 2], "_out_A1FREQ0p001_info0p9.txt")
gwas_MZL <- fread(MZL, header = TRUE, na.strings="NA", sep = " ")
gwas_MZL <- gwas_MZL[,c(19,14,15)]
names(gwas_MZL) <- c("SNP", "MZL.BETA", "MZL.SE")

#merge sumstats

# List of data frames
gwas_list <- list(gwas_CLL, gwas_DLBCL, gwas_FL, gwas_HL, gwas_LPL_WM, gwas_MGUS, gwas_MM, gwas_MZL)

# Perform iterative merging using reduce

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



