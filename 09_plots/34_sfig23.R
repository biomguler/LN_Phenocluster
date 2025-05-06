# ======== Clear environment ========
rm(list = ls())
gc()

# ------------------------
# Load required libraries
# ------------------------
library(data.table)
library(dplyr)

# ------------------------
# Load DrugBank annotation
# ------------------------
load("/omics/groups/OE0136/internal/private/Murat/R_project/LN_annotation/dgi_with_fast_drugbank.Rdata")

drug_info_dt_filter <- drug_info_dt %>%
  filter(!is.na(drugbank_id), atc_code != "") %>%
  distinct() %>%
  as.data.table()
# ------------------------
# Load and annotate DGI data
# ------------------------
dgi_drugs <- fread("dgi_atc.tab", select = c(1,2))  # Must have: drugbank_id, atc_code



# ---------------------
# Define helper: %nin%
# ---------------------
`%nin%` <- Negate(`%in%`)

# ---------------------
# Define ATC codes
# ---------------------
atc_L2_codes <- c("A01", "A04", "A05", "A06", "A07", "A10", "A11", "A14", "A16",
                  "B01", "B02", "B03", "B05",
                  "C01", "C02", "C03", "C04", "C05", "C07", "C08", "C09", "C10",
                  "D01", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11",
                  "G01", "G02", "G03", "G04",
                  "H01", "H02", "H03", "H05",
                  "J01", "J04", "J05", "J07",
                  "L01", "L02", "L03", "L04",
                  "M01", "M02", "M04", "M05",
                  "N01", "N02", "N03", "N04", "N05", "N06", "N07",
                  "P01", "P02", "P03",
                  "R01", "R02", "R03", "R06",
                  "S01", "S02", "S03",
                  "V03", "V04", "V09", "V10")
atc_L1_codes <- sort(unique(substr(atc_L2_codes, 1, 1)))

# ---------------------
# Function to expand ATC codes
# ---------------------
expand_atc <- function(dt) {
  dt_expanded <- dt[, .(atc_code = unlist(strsplit(atc_code, ";"))), by = drugbank_id]
  dt_expanded[, atc_L1 := substr(atc_code, 1, 1)]
  dt_expanded[, atc_L2 := substr(atc_code, 1, 3)]
  return(dt_expanded)
}

# ---------------------
# Load your input data
# ---------------------
# dgi_drugs <- fread("dgi_atc.tab")
# drug_info_dt_filter <- fread("drugbank_background.tab")

# Expand ATC codes
dgi_long <- expand_atc(dgi_drugs)
bg_long <- expand_atc(drug_info_dt_filter)

# ---------------------
# Enrichment function (ORA with two-sided CI)
# ---------------------
run_fisher <- function(code, level = "L2") {
  col <- ifelse(level == "L1", "atc_L1", "atc_L2")
  
  dgi_with_code <- unique(dgi_long[get(col) == code, drugbank_id])
  bg_with_code <- unique(bg_long[get(col) == code, drugbank_id])
  
  bg_all <- unique(bg_long$drugbank_id)
  dgi_all <- unique(dgi_long$drugbank_id)
  bg_only <- setdiff(bg_all, dgi_all)
  
  # Contingency table values
  a <- sum(dgi_all %in% bg_with_code)             # In DGI & has code
  b <- sum(dgi_all %nin% bg_with_code)            # In DGI & no code
  c <- sum(bg_only %in% bg_with_code)             # In BG only & has code
  d <- sum(bg_only %nin% bg_with_code)            # In BG only & no code
  
  if (any(c(a, b, c, d) == 0)) return(NULL)
  
  mat <- matrix(c(a, b, c, d), nrow = 2)
  
  # One-sided p-value for enrichment
  fisher_one <- fisher.test(mat, alternative = "greater")
  
  # Two-sided CI for stable bounds
  fisher_two <- fisher.test(mat, alternative = "two.sided")
  
  data.table(
    ATC_Code = code,
    Level = level,
    A = a, B = b, C = c, D = d,
    OddsRatio = fisher_one$estimate,
    Pvalue = fisher_one$p.value,
    CI_Low = fisher_two$conf.int[1],
    CI_High = fisher_two$conf.int[2]
  )
}

# ---------------------
# Run enrichment for all codes
# ---------------------
results_L1 <- rbindlist(lapply(atc_L1_codes, run_fisher, level = "L1"), fill = TRUE)
results_L2 <- rbindlist(lapply(atc_L2_codes, run_fisher, level = "L2"), fill = TRUE)

# Combine results
results_all <- rbind(results_L1, results_L2, fill = TRUE)

# Adjust p-values
results_all[, FDR := p.adjust(Pvalue, method = "BH")]

# Optional: Keep only enriched terms
results_all <- results_all[OddsRatio > 1]

# Sort by FDR
results_all <- results_all[order(FDR)]



########################################################################################
# ------------------------
# Libraries
# ------------------------
library(data.table)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# ------------------------
# Prepare plot data
# ------------------------
plot_data_all <- copy(results_all)
plot_data_all[, neg_log10_FDR := -log10(FDR)]
plot_data_all[, ATC_L1 := substr(ATC_Code, 1, 1)]

# Define ATC L1 group names
atc_l1_labels <- c(
  A = "ALIMENTARY TRACT AND METABOLISM",
  B = "BLOOD AND BLOOD FORMING ORGANS",
  C = "CARDIOVASCULAR SYSTEM",
  D = "DERMATOLOGICALS",
  G = "GENITO URINARY SYSTEM AND SEX HORMONES",
  H = "SYSTEMIC HORMONAL PREPARATIONS",
  J = "ANTIINFECTIVES FOR SYSTEMIC USE",
  L = "ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS",
  M = "MUSCULO-SKELETAL SYSTEM",
  N = "NERVOUS SYSTEM",
  P = "ANTIPARASITIC PRODUCTS",
  R = "RESPIRATORY SYSTEM",
  S = "SENSORY ORGANS",
  V = "VARIOUS"
)

plot_data_all[, Group := atc_l1_labels[ATC_L1]]

# ------------------------
# Color palette
# ------------------------
atc_groups <- unique(plot_data_all$Group)
color_generator <- colorRampPalette(RColorBrewer::brewer.pal("Paired", "Set1"))
group_palette <- setNames(color_generator(length(atc_groups)), atc_groups)

# ------------------------
# Plot 1: OR vs -log10(FDR)
# ------------------------
p1 <- ggplot(plot_data_all, aes(x = neg_log10_FDR, y = OddsRatio, color = Group)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = group_palette) +
  labs(
    x = expression(-log[10]~"(FDR)"),
    y = "Odds Ratio (OR)",
    color = "ATC Level 1 Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(color = "black", face = "bold")
  )

# ------------------------
# Plot 2: Significant only bar plot
# ------------------------
plot_data_sig <- plot_data_all[FDR < 0.05]
plot_data_sig[, ATC_Code := factor(ATC_Code, levels = plot_data_sig[order(A)]$ATC_Code)]

# Custom labels
custom_labels <- c(
  "L" = "L = ANTINEOPLASTIC-IMMUNOMODULATING",
  "L01" = "L01 = ANTINEOPLASTIC",
  "A" = "A = ALIMENTARY TRACT-METABOLISM",
  "V" = "V = VARIOUS",
  "J" = "J = ANTIINFECTIVES",
  "J01" = "J01 = ANTIBACTERIALS"
)
plot_data_sig[, Label := ifelse(as.character(ATC_Code) %in% names(custom_labels),
                                custom_labels[as.character(ATC_Code)],
                                as.character(ATC_Code))]

p2 <- ggplot(plot_data_sig, aes(x = A, y = ATC_Code, fill = Group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = Label), hjust = -0.05, color = "black", size = 3) +
  scale_fill_manual(values = group_palette) +
  labs(
    x = "Count of DGI Drugs (A)",
    y = "Significant ATC Categories (FDR < 0.05)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(color = "black", face = "bold")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)))

# ------------------------
# Combine plots
# ------------------------
combined_plot <- p1 / p2 + plot_layout(heights = c(1, 1))

# Bubbles Plot
# ------------------------
# Prepare data
# ------------------------
plot_data <- copy(results_all)
plot_data[, neg_log10_FDR := -log10(FDR)]
plot_data[, proportion := A / (A + B)]
plot_data[, ATC_L1 := substr(ATC_Code, 1, 1)]

# Make sure ATC_L1 is a factor with a defined order
l1_levels <- sort(unique(plot_data$ATC_L1))
plot_data[, ATC_L1 := factor(ATC_L1, levels = l1_levels)]

# ------------------------
# Bubble plot
# ------------------------
library(ggplot2)

library(ggrepel)

# Label only significant codes
plot_data_sig <- plot_data[FDR < 0.05]

bubble_plot <- ggplot(plot_data, aes(x = ATC_L1, y = OddsRatio)) +
  geom_jitter(aes(size = proportion, color = neg_log10_FDR), width = 0.25, alpha = 0.9) +
  geom_text_repel(
    data = plot_data_sig,
    aes(label = ATC_Code),
    size = 4,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 1,
    segment.size = 0.2,
    segment.color = "gray50",
    direction = "y"
  ) +
  geom_vline(xintercept = seq(1.5, length(l1_levels) - 0.5, by = 1),
             linetype = "dashed", color = "gray60") +
  scale_color_gradient(low = "skyblue", high = "darkred", name = expression(-log[10]~FDR)) +
  scale_size(range = c(2.5, 10), name = "Proportion (A / A + B)") +
  labs(
    x = "ATC Level 1",
    y = "OR",
    title = "ATC Code Enrichment (L1 & L2, Labeled FDR < 0.05)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

# Show the updated plot
print(bubble_plot)



# ------------------------
# Save results
# ------------------------

ggsave(
  filename = "plot1_ATC_enrichment_plot.tiff",
  plot = combined_plot,
  width = 42,
  height = 29,
  bg="white",
  dpi = 600,
  units = "cm",
  device = "tiff")

ggsave(
  filename = "plot2_ATC_enrichment_plot.tiff",
  plot = bubble_plot,
  width = 12,
  height = 8,
  bg="white",
  dpi = 600,
  units = "in",
  device = "tiff")


fwrite(results_all, "ATC_enrichment_results.tab", sep = "\t")
