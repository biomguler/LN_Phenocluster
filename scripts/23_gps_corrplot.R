#############################################################################
#Title:    GPS-GEV and LDSC results corrplots
#Function: Create custom corrplots
#Author:   Murat Guler
#Time:     May 17th 2024
#############################################################################
rm(list=ls())
gc()
#############################################################################
# Load the required library
library(reshape2)
library(data.table)
library(tidyverse)
library(corrplot)

setwd("/your/input/dir")
# TGPS output formating
system("dos2unix gps.txt")

df <- fread("/your/gps/output/gps.txt")

df <- df %>% rename(gps_scale=`GPS_Min-Max scaling`) 
df$Pheno1 <- sub('LM','LN',df$Pheno1)
df$Pheno2 <- sub('LM','LN',df$Pheno2)

# Create matrix for gps values
gps_matrix <- acast(df, Pheno1 ~ Pheno2, value.var = "gps", mean) # because of permutation Pheno1-Pheno2 and Pheno2-Pheno1 results slightly different and we took mean

ldsc_matrix <- acast(df, Pheno1 ~ Pheno2, value.var = "rg", mean)

# Create matrix for pval values
pval_matrix <- acast(df, Pheno1 ~ Pheno2, value.var = "pval", mean)

ldsc_pval_matrix <- acast(df, Pheno1 ~ Pheno2, value.var = "p", mean)
# Create matrix for gps_scale values
gps_scale_matrix <- acast(df, Pheno1 ~ Pheno2, value.var = "gps_scale", mean)


#First plot: lower triangle
tiff("gps_ldsc.tiff", width = 22, height = 20, units ="cm", res = 600)
corrplot(as.matrix(gps_scale_matrix), is.corr = FALSE, hclust.method = "complete",
         p.mat = pval_matrix, sig.level = 4.76e-4, col = COL2('RdBu', 10),
         insig = "label_sig", pch.cex = 1.2, pch.col = "yellow",
         type = "lower", diag = TRUE, tl.cex=0.8, tl.srt = 60, tl.pos = "lt")

# Second plot: upper triangle
corrplot(as.matrix(ldsc_matrix), is.corr = FALSE, hclust.method = "complete", 
         p.mat = ldsc_pval_matrix, sig.level = 4.76e-4, col = COL2('RdBu', 10),
         insig = "label_sig", pch.cex = 1.2, pch.col = "yellow",
         type = "upper", diag = FALSE, add = TRUE, tl.pos = "n", na.label = "NA", na.label.col = "blue")


dev.off()