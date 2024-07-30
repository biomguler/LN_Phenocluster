###########################################################
# LM cluster analysis
# Murat Guler
# Initiate date: 
# Current date: 
###########################################################
rm(list=ls())
gc()
###########################################################
# Install required packages

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# List of packages to install and load
packages <- c("data.table", "colorspace", "dplyr", "dendextend", "gplots", "tidyverse")

# Install and load each package
for (package in packages) {
  install_if_missing(package)
}

# Load required packages
library(data.table)
library(colorspace)
library(dplyr)
library(dendextend)
library(gplots)
library(tidyverse)

# Approved drugs data
# Upload data file
drug_data <- fread("SData2.txt", sep = "\t")

# Modify format of data
drug_data <- drug_data %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Drug')

drug_data <- drug_data %>% mutate_if(is.numeric,as.logical)

# Create dendrogram

# Dendrogram for drugs
dend_r <-  drug_data %>% dist(method = "binary") %>% hclust(method = "ward.D2") %>% as.dendrogram %>% ladderize

# Dendrogram for LNs
dend_c <- t(drug_data) %>% dist(method = "binary") %>% hclust(method = "ward.D2") %>% as.dendrogram %>% ladderize%>% color_branches(k=3)

# Plot heatmap
tiff("drug_plot.tiff", width = 8, height = 10, units = "in", res = 300)
gplots::heatmap.2(as.matrix(drug_data-1), 
                  main = "Approved Drugs",
                  srtCol = 70,
                  Rowv = dend_r,
                  Colv = dend_c,
                  labRow= TRUE,
                  trace="none",
                  margins =c(6,3),      
                  key.xlab = "no / yes",
                  ylab="Drugs",
                  xlab="LN Subtypes",
                  denscol = "grey",
                  scale="none",
                  density.info = "none",
                  key = FALSE,
                  symkey=FALSE,
                  adjRow =c(0,NA)
)
dev.off()

# Somatic mutation data
# Upload data file
somatic_data <- fread("SData1.txt", sep = "\t")

# Modify format of data
somatic_data <- somatic_data %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Gene')

somatic_data <- somatic_data %>% mutate_if(is.numeric,as.logical)

# Create dendrogram

# Dendrogram for genes
dend_r <-  somatic_data %>% dist(method = "binary") %>% hclust(method = "ward.D2") %>% as.dendrogram %>% ladderize

# Dendrogram for LNs
dend_c <- t(somatic_data) %>% dist(method = "binary") %>% hclust(method = "ward.D2") %>% as.dendrogram %>% ladderize%>% color_branches(k=3)

# Plot heatmap
tiff("somatic_plot.tiff", width = 8, height = 10, units = "in", res = 300)
gplots::heatmap.2(as.matrix(somatic_data-1), 
                  main = "Somatic Mutations",
                  srtCol = 90,
                  Rowv = dend_r,
                  Colv = dend_c,
                  labRow= TRUE,
                  trace="none",
                  margins =c(9,3),      
                  key.xlab = "no / yes",
                  ylab="Mutated Genes",
                  xlab="LN Subtypes",
                  denscol = "grey",
                  scale="none",
                  density.info = "none",
                  na.color = "White",
                  key = FALSE,
                  symkey=FALSE,
                  adjRow =c(0,NA)
)
dev.off()
#End (Line:117)
