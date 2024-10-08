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


PROCESS	/LM_GWAS/MGUS/MGUS_out_A1FREQ0p001_info0p9.txt
PROCESS	/LM_GWAS/MM/MM_out_A1FREQ0p001_info0p9.txt


#for the final meta-analysis of all 8 samples only output results if the


OUTFILE METAL_A1FREQ0p001_MMMGUS_ .tbl
ANALYZE HETEROGENEITY

QUIT
