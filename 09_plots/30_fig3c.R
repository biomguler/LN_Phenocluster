#!/usr/bin/env Rscript
#############################################################################
# Title:    create corrplots ldsc
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
  library(corrplot)
  library(data.table)
})

#############################################################################
# Read the inputs
#############################################################################
df <- fread("ldsc_results.tab")

# Define your column vectors
single_col <- c("CLL", "DLBCL", "FL", "HL", "LPL.WM", "MGUS", "MM", "MCL", "MZL", "PTCL")
phenocluster_col <- c("Cell.B", "Cell.P", "Drug.G1", "MM.MGUS", "Soma.G1", "Soma.G2", "LN")

# Bonferroni correction threshold
bonferroni_threshold <- 0.05 / 10

# Filter single-subtypes data
single_data <- df[p1 %in% single_col & p2 %in% single_col]

# Create empty matrix
mat_single <- matrix(NA, nrow = length(single_col), ncol = length(single_col), dimnames = list(single_col, single_col))
pval_single <- matrix(NA, nrow = length(single_col), ncol = length(single_col), dimnames = list(single_col, single_col))

# Fill matrix
for(i in 1:nrow(single_data)){
  row <- single_data[i]
  mat_single[row$p1, row$p2] <- row$rg
  mat_single[row$p2, row$p1] <- row$rg
  
  pval_single[row$p1, row$p2] <- row$p
  pval_single[row$p2, row$p1] <- row$p
}

diag(mat_single) <- 1

mat_single[mat_single > 1.2 | mat_single < -1.2] <- NA
mat_single[mat_single > 1 & mat_single <= 1.2] <- 1
mat_single[mat_single < -1 & mat_single >= -1.2] <- -1

# Filter phenocluster data (rows single, columns phenocluster)
pheno_data <- df[p1 %in% single_col & p2 %in% phenocluster_col]

mat_pheno <- matrix(NA, nrow = length(single_col), ncol = length(phenocluster_col), dimnames = list(single_col, phenocluster_col))
pval_pheno <- matrix(NA, nrow = length(single_col), ncol = length(phenocluster_col), dimnames = list(single_col, phenocluster_col))

for(i in 1:nrow(pheno_data)){
  row <- pheno_data[i]
  mat_pheno[row$p1, row$p2] <- row$rg
  pval_pheno[row$p1, row$p2] <- row$p
}

# Adjust values greater than abs(1)
mat_pheno[mat_pheno > 1.2 | mat_pheno < -1.2] <- NA
mat_pheno[mat_pheno > 1 & mat_pheno <= 1.2] <- 1
mat_pheno[mat_pheno < -1 & mat_pheno >= -1.2] <- -1

combined_cols <- c(single_col, phenocluster_col)
combined_mat <- matrix(NA, nrow = length(single_col), ncol = length(combined_cols), dimnames = list(single_col, combined_cols))

# Fill lower triangle (single vs single)
combined_mat[single_col, single_col] <- mat_single

# Fill upper triangle (single vs phenocluster)
combined_mat[single_col, phenocluster_col] <- mat_pheno

# Visualization with diagonal solid-filled (blue), no p-values
tiff("corrplot.tiff", width=10, height=6, units="in", res=600, compression="none")

corrplot(combined_mat,
         method = "circle",
         type = "full",
         is.corr = FALSE,
         na.label = "NA",
         tl.col = "black",
         tl.cex = 0.8,
         diag = FALSE,
         mar = c(0,0,1,0),
         tl.srt = 45)

# Add p-values scientifically, highlighting highly significant ones in red
for(i in 1:nrow(combined_mat)){
  for(j in 1:ncol(combined_mat)){
    if(!is.na(combined_mat[i,j])){
      # Identify diagonal cells in single-subtype (first 10 columns)
      is_diagonal_single <- (rownames(combined_mat)[i] == colnames(combined_mat)[j]) && (j <= length(single_col))
      
      if(is_diagonal_single){
        # Solid fill diagonal cells with blue color (no text)
        symbols(j, nrow(combined_mat)-i+1, squares = 0.7, add = TRUE, inches = FALSE,
                fg = "darkblue", bg = "darkblue")
      } else {
        # Select correct p-value matrix
        pvalue <- if(j <= length(single_col)){
          pval_single[rownames(combined_mat)[i], colnames(combined_mat)[j]]
        } else {
          pval_pheno[rownames(combined_mat)[i], colnames(combined_mat)[j]]
        }
        
        if(!is.na(pvalue)){
          # Format p-value scientifically
          p_label <- formatC(pvalue, format = "e", digits = 1)
          
          # Set text color based on significance
          txt_col <- ifelse(pvalue < 0.05/10, "red", "black")
          
          # Add p-value text to plot
          text(j, nrow(combined_mat)-i+1, p_label, col = txt_col, cex = 0.6)
        }
      }
    }
  }
}


dev.off() # close device



