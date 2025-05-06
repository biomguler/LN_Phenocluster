# Load CMplot and gap for lambda GC

library(CMplot)
library(gap)
library(data.table)
library(tidyverse)

# Clear environment
rm(list = ls())
gc()

# Initialize lambda table
lambda_results <- data.frame(trait = character(), lambda_gc = numeric(), stringsAsFactors = FALSE)
# Read sumstats
a <- fread("Combined_GWAS_MultiTrait_for_CMplot.txt")

# Loop over each P-value column (starting from column 4)
for (col in 4:ncol(a)) {
  
  trait <- colnames(a)[col]
  trait_data <- a %>% select("SNP", "CHR", "POS", all_of(trait))
  colnames(trait_data) <- c("SNP", "CHR", "BP", "P")  # standardize for CMplot
  
  # Remove missing p-values
  trait_data <- trait_data[!is.na(trait_data$P), ]
  
  # Skip if no valid p-values
  if (nrow(trait_data) < 100) {
    cat("Skipping", trait, "- not enough data.\n")
    next
  }
  
  # --- λGC calculation ---
  p_clean <- trait_data$P[trait_data$P > 0 & trait_data$P < 1]
  lambda_gc <- gap::gc.lambda(p_clean)
  lambda_results <- rbind(lambda_results, data.frame(trait = trait, lambda_gc = round(lambda_gc, 3)))
  
  # --- QQ Plot ---
  CMplot(trait_data,
         plot.type = "q", 
         box = FALSE, 
         conf.int = TRUE,
         file = "tiff", 
         file.name = paste0("qq", trait),
         dpi = 600, 
         width = 6, 
         height = 6, 
         verbose = TRUE)
  
  # --- Manhattan Plot ---
  CMplot(trait_data,
         type = "p", 
         plot.type = "m", 
         LOG10 = TRUE,
         threshold = 5e-8, 
         chr.labels.angle = 45,
         file = "tiff", 
         file.name = paste0("manhattan", trait),
         dpi = 600, 
         width = 14, 
         height = 6, 
         verbose = TRUE)
}

# ASSET 1-sided plots and lambda
loadGWASData <- function(filepath, phenotype) {
  df <- fread(filepath, header = TRUE, select = c("SNP", "CHR", "POS", "PVAL"))
  df <- df %>%
    dplyr::rename(!!paste0("PVAL_", phenotype) := PVAL)
  return(df)
}

asset_path <- "/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/asset/replication/analysis/sumstats/clean_sumstats/METAL_results/single/asset/merged/rep_asset_1sided_fuma.gz"

asset_sum <- loadGWASData(asset_path, "ASSET")

colnames(asset_sum) <- c("SNP", "CHR", "BP", "P")

p_clean <- asset_sum$P[asset_sum$P > 0 & asset_sum$P < 1]
lambda_gc <- gap::gc.lambda(p_clean)
lambda_results <- rbind(lambda_results, data.frame(trait = "ASSET1-sided", lambda_gc = round(lambda_gc, 3)))

# --- QQ Plot ---
CMplot(asset_sum,
       plot.type = "q", 
       box = FALSE, 
       conf.int = TRUE,
       file = "tiff", 
       file.name = paste0("qq", "Asset"),
       dpi = 600, 
       width = 6, 
       height = 6, 
       verbose = TRUE)

# --- Manhattan Plot ---
CMplot(asset_sum,
       type = "p", 
       plot.type = "m", 
       LOG10 = TRUE,
       threshold = 5e-8, 
       chr.labels.angle = 45,
       file = "tiff", 
       file.name = paste0("manhattan", "Asset"),
       dpi = 600, 
       width = 14, 
       height = 6, 
       verbose = TRUE)



# Save λGC values


write.table(lambda_results, file = "lambda_gc_values.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)




