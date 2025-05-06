#!/usr/bin/env Rscript

# Clear environment
rm(list = ls())
gc()

# Load packages
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggnewscale)
  library(readr)
  library(stringr)
  library(patchwork)
  library(cowplot)
})

# ==== DEFINE FIXED COLORS ====
color_single_col <- c("CLL", "DLBCL", "FL", "HL", "LPL-WM", "MGUS", "MM", "MCL", "MZL", "PTCL")
color_phenocluster_col <- c("Cell-B", "Cell-P", "Drug-G1", "MM-MGUS", "Soma-G1", "Soma-G2", "LN", "ASSET")

color_palette <- c(
  "blue3", "#ff7f0e", "#2ca02c", "#E41A1C", "#9467bd", "#8c564b", "#E7298A", "#7f7f7f", "#17becf", "#bcbd22",
  "#393b79", "#637939", "#843c39", "#31a354", "#3182bd", "#f03b20", "#6a51a3", "firebrick4"
)

named_colors <- setNames(color_palette, c(color_single_col, color_phenocluster_col))

# ==== FUNCTION TO PREP TILEPLOT ====
make_tileplot <- function(input_data, legend = TRUE) {
  
  colnames(input_data)[1] <- "Locus"
  
  # Define column categories
  pheno_cols <- c("Phenocluster", "Subtypes", "Novel")
  evidence_cols <- c("FLAMES", "MAGMA", "OT.L2G", "VEP", "eQTL", "sceQTL", "pQTL", "sQTL", "tuQTL")
  score_col <- "Score"
  all_cols <- c(pheno_cols, evidence_cols, score_col)
  
  # -- PHENOTYPE tiles
  phenotype_tile <- input_data %>%
    select(Locus, all_of(pheno_cols)) %>%
    pivot_longer(-Locus, names_to = "Category", values_to = "Value") %>%
    filter(Value != "") %>%
    mutate(Value = strsplit(Value, ",")) %>%
    mutate(Value = lapply(Value, function(v) sort(trimws(v)))) %>%
    unnest(Value) %>%
    mutate(color = named_colors[Value]) %>%
    group_by(Locus, Category) %>%
    mutate(pheno_count = n(), pheno_index = row_number()) %>%
    ungroup()
  
  gene_levels <- rev(unique(input_data$Locus))
  gene_map <- setNames(seq_along(gene_levels), gene_levels)
  phenotype_tile$y_num <- gene_map[phenotype_tile$Locus]
  
  phenotype_tile <- phenotype_tile %>%
    mutate(
      pheno_width = 0.9 / pheno_count,
      xmin = as.numeric(factor(Category, levels = all_cols)) - 0.45 + (pheno_index - 1) * pheno_width,
      xmax = xmin + pheno_width,
      ymin = y_num - 0.35,
      ymax = y_num + 0.35
    )
  
  # -- EVIDENCE ✓ tiles
  evidence_tile <- input_data %>%
    select(Locus, all_of(evidence_cols)) %>%
    pivot_longer(-Locus, names_to = "Category", values_to = "Value") %>%
    mutate(Status = ifelse(Value == "Yes", "✓", ""),
           y_num = gene_map[Locus])
  
  # -- Score tiles
  scores <- input_data %>%
    select(Locus, Score = !!sym(score_col)) %>%
    mutate(Score = as.numeric(Score), y_num = gene_map[Locus])
  
  # -- Background grid
  grid_lines <- expand.grid(
    Category = all_cols,
    Locus = gene_levels,
    stringsAsFactors = FALSE
  ) %>%
    mutate(y_num = gene_map[Locus])
  
  # === PLOT ===
  p <- ggplot() +
    geom_tile(data = grid_lines,
              aes(x = Category, y = y_num),
              width = 0.9, height = 0.7,
              fill = NA, color = "gray70", linewidth = 0.3) +
    
    # Colored phenotype rectangles
    geom_rect(data = phenotype_tile,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Value),
              color = "white") +
    
    scale_fill_manual(
      name = "Phenotypes",
      values = named_colors,
      guide = if (legend) guide_legend(ncol = 1) else "none"
    ) +
    new_scale_fill() +
    
    # Evidence ✓ tiles
    geom_tile(data = evidence_tile,
              aes(x = Category, y = y_num),
              width = 0.9, height = 0.7,
              fill = "white", color = "gray60") +
    geom_text(data = evidence_tile,
              aes(x = Category, y = y_num, label = Status),
              size = 3, color = "black") +
    
    # Score tiles
    geom_tile(data = scores,
              aes(x = "Score", y = y_num, fill = Score),
              width = 0.9, height = 0.7,
              color = "black") +
    geom_text(data = scores,
              aes(x = "Score", y = y_num, label = round(Score, 2)),
              size = 3, color = "white") +
    
    scale_fill_gradient(name = "Score", low = "lightblue", high = "darkblue", guide = if (legend) "colourbar" else "none") +
    
    scale_y_continuous(breaks = seq_along(gene_levels), labels = gene_levels) +
    scale_x_discrete(limits = all_cols) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 30, face = "bold", hjust = 1, size = 10),
      axis.text.y = element_text(face = "bold", size = 10),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = if (legend) "right" else "none"
    ) +
    coord_fixed(ratio = 0.4)
  
  return(p)
}

# ==== LOAD AND SPLIT ====
input <- read.delim("genes.tsv", stringsAsFactors = FALSE)
part1 <- input[1:68, ]
part2 <- input[69:136, ]

# ==== PLOTS ====
p1 <- make_tileplot(part1, legend = FALSE)
p2 <- make_tileplot(part2, legend = FALSE)
legend_dummy <- make_tileplot(part1, legend = TRUE) + theme(legend.position = "right") + theme_void()
legend_only <- cowplot::get_legend(legend_dummy)

# ==== COMBINE & SAVE ====
final_plot <- cowplot::plot_grid(p1, p2, legend_only, nrow = 1, rel_widths = c(1, 1, 0.3))

ggsave("genes_tileplot_combined_with_evidence.tiff", plot = final_plot, width = 65, height = 40, units = "cm", dpi = 600, bg = "white")
ggsave("genes_tileplot_combined_with_evidence.pdf", plot = final_plot, width = 65, height = 40, units = "cm", dpi = 600, device = cairo_pdf)
