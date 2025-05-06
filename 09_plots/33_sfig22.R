#!/usr/bin/env Rscript
#############################################################################
# Title:    CMplot part2
# Author:   Murat Guler
# Date:     December 3rd, 2024
#############################################################################

# Clear environment
rm(list = ls())
gc()

#############################################################################
# Required packages
#############################################################################
suppressMessages({
  library(CMplot)
  library(data.table)
  library(dplyr)
  library(RColorBrewer)
})

# Define colors
color_single_col <- c("CLL", "DLBCL", "FL", "HL", "LPL-WM", "MGUS", "MM", "MCL", "MZL", "PTCL")
color_phenocluster_col <- c("Cell-B", "Cell-P", "Drug-G1", "MM-MGUS", "Soma-G1", "Soma-G2", "LN")

# Combined color palette (18 colors)
color_palette <- c(
  "blue3",       # CLL
  "#ff7f0e",     # DLBCL
  "#2ca02c",     # FL
  "#E41A1C",     # HL
  "#9467bd",     # LPL-WM
  "#8c564b",     # MGUS
  "#E7298A",     # MM
  "#7f7f7f",     # MCL
  "#17becf",     # MZL
  "#bcbd22",     # PTCL
  "#393b79",     # Cell-B
  "#637939",     # Cell-P
  "#843c39",     # Drug-G1
  "#31a354",     # MM-MGUS
  "#3182bd",     # Soma-G1
  "#f03b20",     # Soma-G2
  "#6a51a3",     # LN
  "firebrick4"   # ASSET
)

# Assign names to the color palette
named_colors <- setNames(color_palette, c(color_single_col, color_phenocluster_col, "ASSET"))



# Trait groups
single_col <- c("CLL", "DLBCL", "FL", "HL", "LPL-WM", "MGUS", "MM", "MCL", "MZL", "PTCL")
phenocluster_col <- c("Cell-B", "Cell-P", "Drug-G1", "MM-MGUS", "Soma-G1", "Soma-G2", "LN")
# Load combined data (assuming this was prepared before)
data <- fread("Combined_GWAS_MultiTrait_for_CMplot.txt")
single_pval_cols <- paste0("PVAL_", single_col)
phenocluster_pval_cols <- paste0("PVAL_", phenocluster_col)

single_data2merge <- data[, c("SNP", "CHR", "POS", single_pval_cols), with = FALSE]
phenocluster_data2merge <- data[, c("SNP", "CHR", "POS", phenocluster_pval_cols), with = FALSE]



loadGWASData <- function(filepath, phenotype) {
  df <- fread(filepath, header = TRUE, select = c("SNP", "CHR", "POS", "PVAL"))
  df <- df %>%
    dplyr::rename(!!paste0("PVAL_", phenotype) := PVAL)
  return(df)
}

asset_path <- "/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/asset/replication/analysis/sumstats/clean_sumstats/METAL_results/single/asset/merged/rep_asset_1sided_fuma.gz"

asset_sum <- loadGWASData(asset_path, "ASSET")

single_combined_data <- merge(asset_sum, single_data2merge, by = c("SNP", "CHR", "POS"), all = TRUE)
phenocluster_combined_data <- merge(asset_sum, phenocluster_data2merge, by = c("SNP", "CHR", "POS"), all = TRUE)
#############################################################################
# Output combined data to file
#############################################################################
single_output_file <- "single_ASSET_Combined_GWAS_MultiTrait_for_CMplot.txt"
phenocluster_output_file <- "phenocluster_ASSET_Combined_GWAS_MultiTrait_for_CMplot.txt"

fwrite(single_combined_data, single_output_file, sep = "\t", na = "NA", quote = FALSE)
fwrite(phenocluster_combined_data, phenocluster_output_file , sep = "\t", na = "NA", quote = FALSE)
message("ASSET_Combined GWAS data saved to: ", c(single_output_file, phenocluster_output_file ))


# Load combined data (assuming this was prepared before)
single_data <- single_combined_data
phenocluster_data <- phenocluster_combined_data

# Define traits for this plot: 1 ASSET + 10 phenoclusters
single_traits <- c("ASSET", single_col)
phenocluster_traits <- c("ASSET", phenocluster_col)

# Build the column names to extract from the data
single_pval_cols <- paste0("PVAL_", single_traits)
phenocluster_pval_cols <- paste0("PVAL_", phenocluster_traits)

 
# Check if all columns exist
if (!all(single_pval_cols %in% colnames(single_data))) {
  warning("Skipping ", single_traits, " because one or more PVAL columns are missing.")
  next
}

if (!all(phenocluster_pval_cols %in% colnames(phenocluster_data))) {
  warning("Skipping ", phenocluster_traits, " because one or more PVAL columns are missing.")
  next
}


# Assign trait names to the colors 
single_signal_colors <- named_colors[single_traits]
phenocluster_signal_colors <- named_colors[phenocluster_traits]
# Optional: Rename columns for prettier legend
setnames(single_data, old = single_pval_cols, new = single_traits)
setnames(phenocluster_data, old = phenocluster_pval_cols, new = phenocluster_traits)  

# Generate ASSET vs single
CMplot(single_data,
       plot.type = "m",
       col = "grey",
       multraits = TRUE,
       threshold = 5e-8,
       ylim = c(0,80),
       threshold.lty = 1,
       threshold.lwd = c(1, 1),
       threshold.col = c("black", "grey"),
       amplify = TRUE,
       chr.den.col = NULL,
       signal.col = single_signal_colors,
       signal.cex = 0.8,
       file = "tiff",
       file.name = paste0("Manhattan_", "ASSET", "_vs_single"),
       dpi = 600,
       file.output = TRUE,
       verbose = TRUE,
       points.alpha = 225,
       legend.ncol = 11,
       legend.pos = "middle")

# Generate ASSET vs Phenocluster
CMplot(phenocluster_data,
       plot.type = "m",
       col = "grey",
       multraits = TRUE,
       threshold = 5e-8,
       ylim = c(0,80),
       threshold.lty = 1,
       threshold.lwd = c(1, 1),
       threshold.col = c("black", "grey"),
       amplify = TRUE,
       chr.den.col = NULL,
       signal.col = phenocluster_signal_colors,
       signal.cex = 0.8,
       file = "tiff",
       file.name = paste0("Manhattan_", "ASSET", "_vs_phenocluster"),
       dpi = 600,
       file.output = TRUE,
       verbose = TRUE,
       points.alpha = 225,
       legend.ncol = 8,
       legend.pos = "middle")

