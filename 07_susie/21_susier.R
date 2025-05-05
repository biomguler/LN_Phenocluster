#!/usr/bin/env Rscript
#############################################################################
# Title:    susie runner
# Author:   Murat Guler
# Time:     December 3rd 2024
#############################################################################
rm(list = ls())
gc()
#############################################################################
# Required packages
#############################################################################
suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(susieR)
  library(arrow)
  library(ggpubr)
  library(RColorBrewer)
})

#############################################################################
# Parse command-line arguments
#############################################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 04_susie.R <finemap_file_path> <gene_file_path>")
}

# Arguments:
#   finemap_file_path: path to finemap_file.

finemap_file <- args[1]
gene_file <- args[2]
message("Using Fine-map file: ", finemap_file)
message("Using Gene file: ", gene_file)
#############################################################################
# Step 1: Load fine-map file
#############################################################################
df_finemap_file <- fread(finemap_file)
if (!file.exists(finemap_file)) {
  stop("Error: Finemap file not found in ", finemap_file)
}

#############################################################################
# Step 2: For each region, subset the data and run susie
#############################################################################

# Check for LSF job array index (if used).
job_index_str <- Sys.getenv("LSB_JOBINDEX")
if (job_index_str != "") {
  region_indices <- as.numeric(job_index_str)
  message("Processing region index: ", region_indices)
} else {
  region_indices <- 1:nrow(df_finemap_file)
}

for (i in region_indices) {
  message("Processing region ", i, " of ", nrow(df_finemap_file))
  
  # Get locus
  locus <- fread(df_finemap_file$sumstatpath[i])
  chr <- as.numeric(df_finemap_file$CHR[i])
  pos <- as.numeric(df_finemap_file$POS[i])
  locus <- locus[locus$CHR == chr & locus$POS < pos + 5e5 & locus$POS > pos - 5e5, ]
  locus <- locus %>% arrange(SNP)
  
  SNP <- data.frame(index = locus$SNP)
  
  # Get LD data
  df_R <- read_feather(df_finemap_file$ldfile_path[i])
  df_R <- as.data.frame(df_R)
  df_R_indexed <- merge(df_R, SNP, by="index")
  rownames(df_R_indexed) <- df_R_indexed[[1]]
  df_R_indexed <- df_R_indexed[, -1] 
  df_R_indexed <- df_R_indexed %>% select(rownames(df_R_indexed))
  
  # Filter locus
  locus_filtered <- locus[locus$SNP %in% rownames(df_R_indexed)]
  
  # Check alignment of sumstats and LD file
  if (!all(locus_filtered$SNP == rownames(df_R_indexed)) | !all(locus_filtered$SNP == colnames(df_R_indexed))) {
    stop("Error: SNPs in sumstats and LD file are not aligned. Please check SNP order.")
  }
  message("Sumstats and LD files aligned.")
  
  # LD matrix
  n_samplesize <- as.numeric(df_finemap_file$n[i])
  ld.mat <- matrix(as.vector(data.matrix(df_R_indexed)), 
                   nrow=nrow(locus_filtered), 
                   ncol=nrow(locus_filtered))
  
  fitted_rss <- susie_rss(bhat = locus_filtered$BETA, 
                          shat = locus_filtered$SE, 
                          R = ld.mat, n = n_samplesize)
  
  locus_filtered$variable <- as.double(rownames(locus_filtered))
  locus_filtered_susie <- merge(locus_filtered, summary(fitted_rss)$vars, by="variable", all.x=TRUE)
  
  lead_snp <- as.character(df_finemap_file$Lead[i])
  lead_snp_variable <- as.double(locus_filtered_susie$variable[locus_filtered_susie$SNP == lead_snp])
  locus_filtered_susie$lead_snp <- as.character(lead_snp)
  
  # Compute R2 with lead SNP
  locus_filtered_susie <- locus_filtered_susie %>% 
    mutate(R2 = (as.numeric(ld.mat[lead_snp_variable, locus_filtered_susie$variable]))^2)
  
  # Get credible sets
  temp.cs <- susie_get_cs(fitted_rss, ld.mat, coverage = 0.95, min_abs_corr = 0.5)
  purity <- as.data.frame(temp.cs$purity)
  coverage_index <- data.frame(cs_index = temp.cs$cs_index, coverage = temp.cs$coverage)
  
  # Get credible set summary
  cs_summary <- summary(fitted_rss)$cs
  
  # If cs_summary is NULL, create an empty placeholder with the same number of rows as purity
  if (is.null(cs_summary)) {
    cs_summary <- data.frame(matrix(NA, ncol = 1, nrow = nrow(purity)))  
    colnames(cs_summary) <- "cs_summary"
  }
  
  # Ensure cs_summary has the same number of rows as purity by adding NA rows at the bottom
  if (nrow(cs_summary) < nrow(purity)) {
    missing_rows <- nrow(purity) - nrow(cs_summary)
    na_df <- as.data.frame(matrix(NA, ncol = ncol(cs_summary), nrow = missing_rows))  # Create NA rows
    colnames(na_df) <- colnames(cs_summary)  # Ensure column names match
    cs_summary <- rbind(cs_summary, na_df)  # Append missing rows safely
  }
  
  # Combine all pieces safely
  cs_info <- cbind(purity, coverage_index, cs_summary)
  
  # Save credible sets & summary
  fwrite(locus_filtered_susie, paste0(df_finemap_file$output_name[i], "_credible_snps.tab"), sep = "\t")
  fwrite(cs_info, paste0(df_finemap_file$output_name[i], "_credible_set_info.tab"), sep = "\t")
  
  # Define color scale for R2
  # Load gene position file
  
  df_genes <- fread(gene_file)
  
  # Filter genes within the locus region
  df_genes_filtered <- df_genes %>%
    filter(CHR == chr & start <= max(locus_filtered_susie$POS) & end >= min(locus_filtered_susie$POS))
  
  # Define color scale for R2
  color_scale <- scale_color_gradient(low = "blue", high = "red")
  
  # Define shape scale for credible sets
  shape_values <- c(16, 17, 18, 15, 3, 8)
  unique_cs <- unique(na.omit(locus_filtered_susie$cs)) 
  shape_map <- setNames(shape_values[seq_along(unique_cs)], unique_cs)
  
  # Get CS shape for lead SNP
  lead_snp_cs <- locus_filtered_susie$cs[locus_filtered_susie$SNP == lead_snp]
  lead_snp_shape <- ifelse(!is.na(lead_snp_cs), shape_map[as.character(lead_snp_cs)], 16)
  
  # Get lead SNP position
  lead_snp_pos <- locus_filtered_susie$POS[locus_filtered_susie$SNP == lead_snp]
  
  # Plot -log10 P-value (Keeps the Legend)
  a <- ggplot(locus_filtered_susie, aes(x = POS, y = -log10(PVAL))) +
    geom_point(aes(color = R2, shape = as.factor(cs)), alpha = 0.8, size = 2) +
    geom_point(data = subset(locus_filtered_susie, SNP == lead_snp), 
               aes(shape = as.factor(cs)), color = "purple", size = 5) +
    geom_label_repel(data = subset(locus_filtered_susie, cs != -1), aes(label = SNP), 
                     size = 4, max.overlaps = 15) +
    geom_vline(xintercept = lead_snp_pos, linetype = "dashed", color = "black", size = 1) +  # Vertical line
    color_scale +
    scale_shape_manual(values = shape_map) +
    labs(title = "Fine-Mapping Regional Plot", 
         y = "-log10(P-value)", 
         color = "R²", 
         shape = "Credible Set") +
    theme_minimal() +
    theme(legend.position = "right",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # Plot Posterior Inclusion Probability (PIP) (Removes Legend & X-Axis)
  b <- ggplot(locus_filtered_susie, aes(x = POS, y = variable_prob)) +
    geom_point(aes(color = R2, shape = as.factor(cs)), alpha = 0.8, size = 2) +
    geom_point(data = subset(locus_filtered_susie, SNP == lead_snp), 
               aes(shape = as.factor(cs)), color = "purple", size = 5) +
    geom_label_repel(data = subset(locus_filtered_susie, cs != -1), aes(label = SNP), 
                     size = 4, max.overlaps = 15) +
    geom_vline(xintercept = lead_snp_pos, linetype = "dashed", color = "black", size = 1) +  # Vertical line
    color_scale +
    scale_shape_manual(values = shape_map) +
    labs(title = "Posterior Inclusion Probability (PIP)", 
         y = "PIP") +  # Remove x-axis label
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())  # Removes x-axis from plot B
  
  # Generate distinct colors using the "Set1" palette from RColorBrewer
  set.seed(123)  # Ensures reproducibility
  num_genes <- length(unique(df_genes_filtered$name))
  
  # Extend the "Set1" palette if there are more than 9 genes
  gene_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(num_genes), unique(df_genes_filtered$name))
  
  # Ensure gene names are properly ordered in y-axis
  df_genes_filtered$name <- factor(df_genes_filtered$name, levels = rev(names(gene_colors)))
  
  # Define min and max position limits to prevent trimming issues
  plot_min <- min(locus_filtered_susie$POS)
  plot_max <- max(locus_filtered_susie$POS)
  
  # Trim gene start/end positions if they exceed the plot range
  df_genes_filtered$start_trimmed <- pmax(df_genes_filtered$start, plot_min)  # Ensures start is within range
  df_genes_filtered$end_trimmed <- pmin(df_genes_filtered$end, plot_max)  # Ensures end is within range
  
  # Plot Gene Positions as Horizontal Lines in a New Panel (C)
  c <- ggplot(df_genes_filtered, aes(xmin = start_trimmed, xmax = end_trimmed, y = name)) +
    geom_segment(aes(x = start_trimmed, xend = end_trimmed, y = name, yend = name, color = name), size = 1) +  # Apply colors to bars
    geom_vline(xintercept = lead_snp_pos, linetype = "dashed", color = "black", size = 1) +  # Vertical line
    scale_x_continuous(limits = c(plot_min, plot_max)) +  # Align x-axis
    scale_color_manual(values = gene_colors, guide = "none") +  # Apply colors & remove legend
    labs(title = "Gene Annotations", x = "Genomic Position", y = "Genes") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8, color = gene_colors[df_genes_filtered$name]))  # Fix color issue
  
  # Save the plot
  output_name <- paste0(df_finemap_file$output_name[i], "_fm_regionalplot.jpeg")
  final_plot <- ggarrange(a, b, c, ncol = 1, nrow = 3, align = "v", heights = c(2, 2, 1))  # Adjust height ratio
  
  ggsave(output_name, plot = final_plot, dpi = 600, width = 10, height = 12, units = "in")
  
  
}
