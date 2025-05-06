#!/usr/bin/env Rscript
#############################################################################
# Title:    create upset and venn plots
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
  library(UpSetR)
  library(VennDiagram)
  library(ggplot2)
  library(data.table)
})

#############################################################################
# Read the inputs
#############################################################################

all <- read.csv("binary_matrix.csv", header=T, sep="," )
groups <- read.csv("summary_binary.csv", header=T, sep=",")


#############################################################################
# Create UpSet Plot 1
#############################################################################

tiff("upset_plot1.tiff", width=13, height=10, units="in", res=600, compression="none") 

upset(all, 
      sets = colnames(all)[-1],
      sets.bar.color = "darkblue",
      matrix.color = "purple",
      main.bar.color = "darkblue",
      order.by = c("freq"),
      show.numbers ="yes",
      mainbar.y.label = "Co-occurrence of independent loci",
      sets.x.label = "Total number of loci",
      text.scale = 2, 
      #group.by ="degree",
      mainbar.y.max = 9, 
      shade.color = "lightgray",
      mb.ratio = c(0.55, 0.45),
      shade.alpha = 0.75,
      matrix.dot.alpha = 0.75, point.size = 2.5,
      set_size.show = T,
      set_size.numbers_size= 10,
      set_size.scale_max=65)

dev.off()

#############################################################################
# Create UpSet Plot 2
#############################################################################

tiff("upset_plot2.tiff", width = 12, height = 8, units = "in", res = 600, compression = "none")

upset(groups, 
      sets = colnames(groups)[-1],
      sets.bar.color = "darkblue",
      matrix.color = "purple",
      main.bar.color = "darkblue",
      order.by = c("freq", "degree"),
      show.numbers = "yes",
      mainbar.y.label = "Co-occurrence of independent loci",
      sets.x.label = "Total number of loci",
      text.scale = 2,
      mb.ratio = c(0.55, 0.45),
      group.by = "degree",
      mainbar.y.max = 30,
      shade.color = "lightgray",
      shade.alpha = 0.75,
      matrix.dot.alpha = 0.75,
      point.size = 2.5,
      set_size.show = TRUE,
      set_size.scale_max=65, 
      set_size.numbers_size= 10)

dev.off()

#############################################################################
# Create Venn Diagram
#############################################################################

# Extract sets
asset_set <- groups$Index[groups$ASSET == 1]
phenocluster_set <- groups$Index[groups$Phenocluster == 1]
single_set <- groups$Index[groups$Single == 1]

# Save Venn Diagram at 600 DPI
tiff("venn_diagram.tiff", width=12, height=12, units="in", res=600, compression="none") 

venn.plot <- draw.triple.venn(
  area1 = length(asset_set),
  area2 = length(phenocluster_set),
  area3 = length(single_set),
  n12 = length(intersect(asset_set, phenocluster_set)),
  n23 = length(intersect(phenocluster_set, single_set)),
  n13 = length(intersect(asset_set, single_set)),
  n123 = length(intersect(intersect(asset_set, phenocluster_set), single_set)),
  category = c("ASSET", "Phenocluster", "Single"),
  fill = c("#fbe5d6", "#c5e0b4", "#a0bcd4"),
  fontfamily = rep("arial", 7),
  label.col = rep("black", 7), 
  cex = rep(2, 7), 
  fontface = rep("plain", 7),
  alpha = rep(0.5, 3),
  cat.cex = 1.5,
  lwd = 2
)

grid.draw(venn.plot)
dev.off()
