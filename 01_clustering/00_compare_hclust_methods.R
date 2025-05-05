#############################################################################
# Title:    Compare hclust algorithms
# Function: Create correlation plots and compute Fowlkes-Mallows Index
# Author:   Murat Guler (murat.guler@dkfz.de, muratgmbg@gmail.com)
# Date:     Jan 23th 2022
# Note:     Please let me know if you have any trouble
#############################################################################
# Free R memory and remove prior environment
rm(list=ls())
gc()

# Increase stack size to handle large dendrograms
options(expressions = 500000)
#############################################################################
# Packages

# Install required packages, if not installed

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# List of packages to install and load
packages <- c("data.table", "colorspace", "dplyr", "dendextend", "gplots", "tidyverse", "magrittr", "parallel")

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
library(parallel)
library(magrittr)

#############################################################################

# Approved drugs data
# Upload data file
drug_data <- fread("SData2.txt", sep = "\t")

# Modify format of data
drug_data <- drug_data %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Drug')

drug_data <- drug_data %>% mutate_if(is.numeric,as.logical)

# Compare hclust algorithms 

hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", 
                    "median", "centroid", "ward.D2")
drug_dendlist <- dendlist()

for(i in seq_along(hclust_methods)) {
  tmp_dend <-  drug_data %>% dist(method = "binary") %>% 
    hclust(method = hclust_methods[i]) %>% as.dendrogram 
  drug_dendlist <- dendlist(drug_dendlist, tmp_dend)
}

names(drug_dendlist) <- hclust_methods

# cophenetic correlation
cophenetic_cors <- cor.dendlist(drug_dendlist)


## Fowlkes-Mallows Index.

FM_cors_k3 <- cor.dendlist(drug_dendlist, method = "FM_index", k = 3)
FM_cors_k4 <- cor.dendlist(drug_dendlist, method = "FM_index", k = 4)
FM_cors_k5 <- cor.dendlist(drug_dendlist, method = "FM_index", k = 5)


# Create corrplots and save

tiff("drugs_correlation_plots.tiff", width = 2400, height = 3000, res = 300)

par(mfrow = c(2, 2))  # 2x2 layout

# Plot correlation 
corrplot::corrplot(cophenetic_cors, "circle", "lower", addCoef.col ='red', number.cex = 0.6)
mtext("cophenetic correlation", side = 4, line = 1, adj = 0.5, cex = 1)

corrplot::corrplot(FM_cors_k3, "circle", "lower", addCoef.col ='red', number.cex = 0.6)
mtext("k = 3", side = 4, line = 1, adj = 0.5, cex = 1)

corrplot::corrplot(FM_cors_k4, "circle", "lower", addCoef.col ='red', number.cex = 0.6)
mtext("k = 4", side = 4, line = 1, adj = 0.5, cex = 1)

corrplot::corrplot(FM_cors_k5, "circle", "lower", addCoef.col ='red', number.cex = 0.6)
mtext("k = 5", side = 4, line = 1, adj = 0.5, cex = 1)


dev.off()

#############################################################################


# Somatic mutation data

# Upload data file
somatic_data <- fread("SData1.txt", sep = "\t")
# Filter number of genes because of data size
somatic_data <- somatic_data %>%
  mutate(Mean = rowMeans(select(., 2:19)))
somatic_data <- somatic_data %>% filter(Mean >= 0.2)
somatic_data  <- somatic_data[, c(1:19)]
# Modify format of data
somatic_data <- somatic_data %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Gene') 
  
somatic_data <- somatic_data %>% mutate_if(is.numeric,as.logical)

# Compare hclust algorithms 

hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", 
                    "median", "centroid", "ward.D2")

somatic_dendlist <- dendlist()

for(i in seq_along(hclust_methods)) {
  tmp_dend <-  somatic_data %>% dist(method = "binary") %>% 
    hclust(method = hclust_methods[i]) %>% as.dendrogram 
  somatic_dendlist <- dendlist(somatic_dendlist, tmp_dend)
}

names(somatic_dendlist) <- hclust_methods

# cophenetic correlation
cophenetic_cors <- cor.dendlist(somatic_dendlist)


## Fowlkes-Mallows Index.

FM_cors_k3 <- cor.dendlist(somatic_dendlist, method = "FM_index", k = 3)
FM_cors_k4 <- cor.dendlist(somatic_dendlist, method = "FM_index", k = 4)
FM_cors_k5 <- cor.dendlist(somatic_dendlist, method = "FM_index", k = 5)


# Create corrplots and save

tiff("somatic_correlation_plots.tiff", width = 2400, height = 3000, res = 300)

par(mfrow = c(2, 2))  # 2x2 layout

# Plot correlation 
corrplot::corrplot(cophenetic_cors, "circle", "lower", addCoef.col ='red', number.cex = 0.6)
mtext("cophenetic correlation", side = 4, line = 1, adj = 0.5, cex = 1)

corrplot::corrplot(FM_cors_k3, "circle", "lower", addCoef.col ='red', number.cex = 0.6)
mtext("k = 3", side = 4, line = 1, adj = 0.5, cex = 1.5)

corrplot::corrplot(FM_cors_k4, "circle", "lower", addCoef.col ='red', number.cex = 0.6)
mtext("k = 4", side = 4, line = 1, adj = 0.5, cex = 1.5)

corrplot::corrplot(FM_cors_k5, "circle", "lower", addCoef.col ='red', number.cex = 0.6)
mtext("k = 5", side = 4, line = 1, adj = 0.5, cex = 1.5)


dev.off()
#End (Line:170)
