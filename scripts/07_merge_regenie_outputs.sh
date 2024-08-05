#!/bin/bash

# List of regenie phenotypes
phenotypes="Cell_B Cell_P CLL DLBCL Drug_G1 FL HL LM LPL_WM MGUS MM_MGUS MM MZL Soma_G1 Soma_G2"

# Create all chr merged regenie output files for each phenotype
for pheno in $phenotypes; do 
  # Extract the header from the first file
  zcat ukb_step2_LM_chr1_${pheno}.regenie.gz | head -1 > assoc.${pheno}.all.regenie
  
  # Find all matching files and concatenate their contents, excluding the header lines
  for file in *${pheno}.gz; do
    zcat "$file" | tail -n +2
  done >> assoc.${pheno}.all.regenie
done
