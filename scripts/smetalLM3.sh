#Setup
SEPARATOR WHITESPACE
SCHEME STDERR
STDERR SE 
GENOMICCONTROL ON
#THIS SCRIPT EXECUTES AN ANALYSIS OF EIGHT STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Cell_Bp0.9.txt...

#LOAD THE INPUT FILES

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER SNP
ALLELE ALLELE1 ALLELE0
EFFECT BETA
PVALUE P 
FREQ A1FREQ
STDERR SE


PROCESS	/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/CLL/CLL_out_A1FREQ0p001_info0p9.txt
PROCESS	/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/DLBCL/DLBCL_out_A1FREQ0p001_info0p9.txt
PROCESS	/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/FL/FL_out_A1FREQ0p001_info0p9.txt
PROCESS	/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/HL/HL_out_A1FREQ0p001_info0p9.txt
PROCESS	/omics/groups/OE0136/internal/private/Murat/UKB/669373/Genomics/Genotypes/Genotype_Results/Genotype_calls/LM_GWAS/MZL/MZL_out_A1FREQ0p001_info0p9.txt

#for the final meta-analysis of all 8 samples only output results if the

OUTFILE METAL_A1FREQ0p001_CellB_ .tbl
ANALYZE HETEROGENEITY

QUIT
