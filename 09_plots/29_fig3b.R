#!/usr/bin/env Rscript
#############################################################################
# Title:    Tileplot with colored phenotypes and HyPrColocPP score
#############################################################################

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
})

# ==== DEFINE COLORS ====
color_single_col <- c("CLL", "DLBCL", "FL", "HL", "LPL-WM", "MGUS", "MM", "MCL", "MZL", "PTCL")
color_phenocluster_col <- c("Cell-B", "Cell-P", "Drug-G1", "MM-MGUS", "Soma-G1", "Soma-G2", "LN", "ASSET")

color_palette <- c(
  "blue3", "#ff7f0e", "#2ca02c", "#E41A1C", "#9467bd", "#8c564b", "#E7298A", "#7f7f7f", "#17becf", "#bcbd22",
  "#393b79", "#637939", "#843c39", "#31a354", "#3182bd", "#f03b20", "#6a51a3", "firebrick4"
)
named_colors <- setNames(color_palette, c(color_single_col, color_phenocluster_col))

# Read data
input <- read.delim("pleiotropy_assesment.tsv", stringsAsFactors = FALSE)

# Columns to be visualized like phenotype tiles (color)
pheno_cols <- c("Multi.trait.signal",
                "Single.trait.GWS", "Single.trait.GWSu", "Single.reported", "Multi.trait.reported",
                "HyPrColoc", "Consensus.Primary", "Consensus.Supportive", "Overall.subtype.s.")

# Text-only columns
novel_cols <- c("Novel.single", "Novel.multi.trait")

# Score and checkmark columns
score_col <- "HyPrColocPP"
final_col <- "Pleiotropy"

# ==== PHENOTYPE TILE PREP ====
phenotype_tile <- input %>%
  select(Locus, all_of(pheno_cols)) %>%
  pivot_longer(-Locus, names_to = "Category", values_to = "Value") %>%
  filter(Value != "") %>%
  mutate(Value = strsplit(Value, ",")) %>%
  mutate(Value = lapply(Value, function(v) sort(trimws(v)))) %>%
  unnest(Value) %>%
  mutate(
    color = named_colors[Value],
    Category = as.character(Category)
  ) %>%
  group_by(Locus, Category) %>%
  mutate(
    pheno_count = n(),
    pheno_index = row_number()
  ) %>%
  ungroup()

# Map gene levels
gene_levels <- rev(unique(input$Locus))
gene_map <- setNames(seq_along(gene_levels), gene_levels)
phenotype_tile$y_num <- gene_map[phenotype_tile$Locus]

# Positioning
phenotype_tile <- phenotype_tile %>%
  mutate(
    pheno_width = 0.9 / pheno_count,
    xmin = as.numeric(factor(Category, levels = c(pheno_cols, novel_cols, "HyPrColocPP", "Pleiotropy"))) - 0.45 + (pheno_index - 1) * pheno_width,
    xmax = xmin + pheno_width,
    ymin = y_num - 0.35,
    ymax = y_num + 0.35
  )

# ==== NOVEL TEXT TILE ====
novel_text <- input %>%
  select(Locus, all_of(novel_cols)) %>%
  pivot_longer(-Locus, names_to = "Category", values_to = "Text") %>%
  filter(Text != "") %>%
  mutate(y_num = gene_map[Locus])

# ==== SCORE TILE ====
scores <- input %>%
  select(Locus, Score = !!sym(score_col)) %>%
  filter(Score != "") %>%
  mutate(Score = as.numeric(Score), y_num = gene_map[Locus])

# ==== PLEIOTROPY TILE ====
pleio_data <- input %>%
  select(Locus, Pleiotropy) %>%
  mutate(
    Status = ifelse(Pleiotropy == "Yes", "✓", ""),
    y_num = gene_map[Locus]
  )

# ==== GRID BACKGROUND (BORDERS) ====
grid_lines <- expand.grid(
  Category = c(pheno_cols, novel_cols, "HyPrColocPP", "Pleiotropy"),
  Locus = gene_levels,
  stringsAsFactors = FALSE
) %>%
  mutate(y_num = gene_map[Locus])

# ==== FINAL PLOT ====
p <- ggplot() +
  # Grid outlines
  geom_tile(data = grid_lines,
            aes(x = Category, y = y_num),
            width = 0.9, height = 0.7,
            fill = NA, color = "gray70", linewidth = 0.3) +
  
  # Phenotype tiles
  geom_rect(data = phenotype_tile,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Value),
            color = "white") +
  
  scale_fill_manual(
    name = "Phenotypes",
    values = named_colors,
    breaks = names(named_colors),
    guide = guide_legend(ncol = 1)
  ) +
  
  new_scale_fill() +
  
  # Novel columns as text-only tiles
  geom_tile(data = novel_text,
            aes(x = Category, y = y_num),
            width = 0.9, height = 0.7,
            fill = "white", color = "gray60") +
  geom_text(data = novel_text,
            aes(x = Category, y = y_num, label = Text),
            size = 2, color = "black") +
  
  # Score tiles
  geom_tile(data = scores,
            aes(x = "HyPrColocPP", y = y_num, fill = Score),
            width = 0.9, height = 0.7,
            color = "black") +
  geom_text(data = scores,
            aes(x = "HyPrColocPP", y = y_num, label = round(Score, 2)),
            size = 3, color = "white") +
  
  scale_fill_gradient(name = "HyPrColocPP", low = "lightblue", high = "darkblue") +
  
  # Pleiotropy ✓
  geom_tile(data = pleio_data,
            aes(x = "Pleiotropy", y = y_num),
            width = 0.9, height = 0.7,
            fill = "white", color = "gray80") +
  geom_text(data = pleio_data,
            aes(x = "Pleiotropy", y = y_num, label = Status), size = 4) +
  
  # Axis settings
  scale_y_continuous(breaks = seq_along(gene_levels), labels = gene_levels) +
  scale_x_discrete(limits = c(pheno_cols, novel_cols, "HyPrColocPP", "Pleiotropy")) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, face = "bold", hjust = 1, size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 0.4)

ggsave("fig_3c_tileplot.tiff", plot = p, device = "tiff", width = 30, height = 38, units = "cm", dpi = 600, bg = "white")


ggsave("fig_3c_tileplot.pdf", plot = p,
       width = 30, height = 38, units = "cm", dpi = 600,
       device = cairo_pdf)