#!/bin/bash
#this script creates combine bed, bim, fam files  for UK Biobank genotype data which are chr separeted

#main code
#step1: combine all chromosome
#remove any file names as output before
rm -f list_beds.txt

#create a txt file to list all bed, bim fam files each line for per chromosome
#please change XXX to you own UKBB application ID number
for chr in {2..22}; do echo "ukbXXX_c${chr}_b0_v2.bed ukb_snp_chr${chr}_v2.bim ukbXXX_c10_b0_v2_s488156.fam" >> list_beds.txt; done
#note all .fam files are same, so it does not matter 

#use plink1.9 to combine
plink \
  --bed ukbXXX_c1_b0_v2.bed \
  --bim ukb_snp_chr1_v2.bim \
  --fam ukbXXX_c1_b0_v2_s488156.fam \
  --merge-list list_beds.txt \
  --make-bed --out ukb_cal_allChrs

#step2: Check genotype data by using plink2 and obtion list of SNP that are passed QC 
plink2 
  --bfile ukb_cal_allChrs \
  --maf 0.01 \
  --mac 5 \
  --geno 0.1 \
  --hwe 1e-8 \
  --mind 0.1 \
  --write-snplist \
  --autosome \
  --write-samples \
  --no-id-header \
  --out qc_pass
