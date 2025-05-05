#############################################################################
#Title:    Fuma results to functional annotations
#Function: Create summary table
#Author:   Murat Guler
#Time:     May 17th 2024
#############################################################################
#############################################################################
rm(list=ls())
gc()
#############################################################################
#required packages
library(data.table)
library(tidyverse)

# Location of result files
setwd("/your/FUMA_results/dir")

# read ci files

ci_files <- list.files(path = "/your/FUMA_results/dir",
                          pattern = "^ci\\.txt$",
                          full.names = TRUE,
                          include.dirs = TRUE,
                          recursive = TRUE)

ci_names <- basename(dirname(ci_files))

# read magma files

magma_files <- list.files(path = "/your/FUMA_results/dir",
                          pattern = "magma.genes.out$",
                          full.names = TRUE,
                          include.dirs = TRUE,
                          recursive = TRUE)

magma_names <- basename(dirname(magma_files))

# read eqtl files
eqtl_files <- list.files(path = "/your/FUMA_results/dir",
                    pattern = "eqtl.txt$",
                    full.names = TRUE,
                    include.dirs = TRUE,
                    recursive = TRUE)

eqtl_names <- basename(dirname(eqtl_files))

# read ld files
ld_files <- list.files(path = "/your/FUMA_results/dir",
                         pattern = "ld.txt$",
                         full.names = TRUE,
                         include.dirs = TRUE,
                         recursive = TRUE)

ld_names <- basename(dirname(ld_files))

# read ld files
ld_files <- list.files(path = "/your/FUMA_results/dir",
                       pattern = "ld.txt$",
                       full.names = TRUE,
                       include.dirs = TRUE,
                       recursive = TRUE)

ld_names <- basename(dirname(ld_files))

# read snp files
snp_files <- list.files(path = "/your/FUMA_results/dir",
                       pattern = "snps.txt$",
                       full.names = TRUE,
                       include.dirs = TRUE,
                       recursive = TRUE)

snp_names <- basename(dirname(snp_files))



# annotations
annot_files <- list.files(path = "/your/FUMA_results/dir",
                        pattern = "annot.txt$",
                        full.names = TRUE,
                        include.dirs = TRUE,
                        recursive = TRUE)

annot_names <- basename(dirname(annot_files))


# annovar
annov_files <- list.files(path = "/your/FUMA_results/dir",
                          pattern = "annov.txt$",
                          full.names = TRUE,
                          include.dirs = TRUE,
                          recursive = TRUE)

annov_names <- basename(dirname(annov_files))

# genes
genes_files <- list.files(path = "/your/FUMA_results/dir",
                          pattern = "genes.txt$",
                          full.names = TRUE,
                          include.dirs = TRUE,
                          recursive = TRUE)

genes_names <- basename(dirname(genes_files))



################################################################################
# store eqtl 
eqtl <- list()
ld <- list()
snp <- list()
annot <- list()
annov <- list()
genes <- list()
magma <- list()
ci <- list()
# Loop through file
for (i in 1:length(ci_files)) {
  # Read the file
  df <- fread(ci_files[i])
  df<- df %>% mutate(pheno=ci_names[i])
  ci[[i]] <- df
}

for (i in 1:length(magma_files)) {
  # Read the file
  df <- fread(magma_files[i])
  df<- df %>% mutate(pheno=magma_names[i])
  df$FDR <- p.adjust(df$P, method= "fdr")
  magma[[i]] <- df
}

for (i in 1:length(eqtl_files)) {
  # Read the file
  df <- fread(eqtl_files[i])
  df<- df %>% mutate(pheno=eqtl_names[i])
  eqtl[[i]] <- df
  }

for (i in 1:length(ld_files)) {
  # Read the file
  df <- fread(ld_files[i])
  df<- df %>% mutate(pheno=ld_names[i])
  ld[[i]] <- df
}

for (i in 1:length(snp_files)) {
  # Read the file
  df <- fread(snp_files[i])
  df<- df %>% mutate(pheno=snp_names[i])
  snp[[i]] <- df
}

for (i in 1:length(annot_files)) {
  # Read the file
  df <- fread(annot_files[i])
  df<- df %>% mutate(pheno=annot_names[i])
  annot[[i]] <- df
}
for (i in 1:length(annov_files)) {
  # Read the file
  df <- fread(annov_files[i])
  df<- df %>% mutate(pheno=annov_names[i])
  annov[[i]] <- df
}
for (i in 1:length(genes_files)) {
  # Read the file
  df <- fread(genes_files[i])
  df<- df %>% mutate(pheno=genes_names[i])
  genes[[i]] <- df
}

names(eqtl) <- eqtl_names
names(ld) <- ld_names
names(snp) <- snp_names
names(annot) <- annot_names
names(annov) <- annov_names
names(genes) <- genes_names
names(magma) <- magma_names
names(ci) <- ci_names


combined_magma <- rbindlist(magma)
combined_eqtl <- rbindlist(eqtl)
combined_ld <- rbindlist(ld)
combined_snp <- rbindlist(snp, fill=TRUE)
combined_annot <- rbindlist(annot, fill=TRUE)
combined_annov <- rbindlist(annov)
combined_genes <- rbindlist(genes, fill=TRUE)
combined_ci <- rbindlist(ci)

# Filter eqtls

filtered_eqtl <- combined_eqtl %>%
  filter(db %in% c("eQTLGen", "GTEx/v8")) %>% 
  filter(tissue  %in% c("eQTLGen_cis_eQTLs", "Cells_EBV-transformed_lymphocytes", 
                        "Whole_Blood"))
# merge

merged_snp_eqlt <- merge(combined_snp, filtered_eqtl, by = c("uniqID", "pheno"), all = TRUE)


merged_snp_eqlt <- merged_snp_eqlt %>% 
  mutate(ID = paste0(chr.x, ":", pos.x, ":", non_effect_allele, ":", effect_allele))

# arrange
merged_snp_eqlt <- merged_snp_eqlt %>% filter(r2 >= 0.6)

# genes
library(tidyr)
genes_separated <- combined_genes %>%
  separate_rows(IndSigSNPs, sep = ";") %>% rename(IndSigSNP = IndSigSNPs)

# magma filter FDR

magma_filtered_FDR <- combined_magma %>% filter(FDR <= 0.05)

# ci filter tidy
ci_separated <- combined_ci %>%
  separate_rows(genes, sep = ";") %>% filter (`tissue/cell` != "Spleen")

ci_separated <- ci_separated %>%
  separate_rows(genes, sep = ":") %>% filter (!is.na(genes))

library(gprofiler2) # to convert ensg to gene symbol
ci_separated$symbols <- gconvert(ci_separated$genes,organism="hsapiens",target="HGNC",filter_na = F)$target


ci_separated <- ci_separated %>% separate_rows(SNPs, sep = ";")



# write
fwrite(merged_snp_eqlt, "merged_snp_eqlt.txt", sep = "\t")
fwrite(genes_separated, "genemapping_byld.txt", sep= "\t")
fwrite(magma_filtered_FDR, "magma_fdr.txt", sep= "\t")
fwrite(ci_separated, "ci_separated.txt", sep= "\t")


