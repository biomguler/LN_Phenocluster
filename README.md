LN Phenocluster Manuscript Pipeline
-----------------------------

**Project Name:** Lymphoid neoplasms (LN) phenocluster

**Description:** This repository contains the scripts and pipeline used in the LN phenocluster manuscript. Please follow the scripts listed below to generate the analysis results described in the manuscript.

**Citation:** "M. Guler & F. Canzian, Clustering of Lymphoid Neoplasms by Cell of Origin, Somatic Mutation and Drug Usage Profiles: A Multi-trait Genome-Wide Association Study, ***submitted(doi:will be updated)***, 2025"

* * * * *

### Required Softwares

| **Software** | **Version** | **Link/Citation** |
| --- | --- | --- |
| regenie | 3.2.1 | [regenie](https://github.com/rgcgithub/regenie) |
| plink | 1.90b6.21 | [plink](https://www.cog-genomics.org/plink/) |
| plink2 | 2.00a3LM | [plink2](https://www.cog-genomics.org/plink/2.0/) |
| R | 4.3.0 or above | [R](https://www.r-project.org/) |
| METAL | 2020-05-05 | [METAL](https://github.com/statgen/METAL) |
| LDSC | v1.0.1 | [LDSC](https://github.com/bulik/ldsc) |
| GWASLab | v3.4.31 | [GWASLab](https://github.com/Cloufield/gwaslab) |
| FLAMES | v1.1.2 | [FLAMES](https://github.com/Marijn-Schipper/FLAMES) |

### Required R Packages

| **R Package** | **Version** | **Link/Citation** |
| --- | --- | --- |
| ASSET | 2.18.0 | [ASSET](https://bioconductor.org/packages/release/bioc/html/ASSET.html) |
| tidyverse | 2.0.0 | [tidyverse](https://www.tidyverse.org/) |
| ggplot2 | 3.1.3.1 | [gplots](https://github.com/talgalili/gplots) |
| dendextend | 1.17.1 | [dendextend](https://academic.oup.com/bioinformatics/article/31/22/3718/240978) |
| data.table | 1.14.10 | [data.table](https://github.com/Rdatatable/data.table) |
| parallel | 4.3.0 | [parallel](https://www.R-project.org/) |
| magrittr | 2.0.3 | [magrittr](https://CRAN.R-project.org/package=magrittr) |
| colorspace | 2.1.0 | [colorspace](https://doi.org/10.1016/j.csda.2008.11.033) |
| susieR | 0.12.35 | [susieR](https://github.com/stephenslab/susieR) |
| HyPrColoc | 0.0.2 | [HyPrColoc](https://github.com/jrs95/hyprcoloc) |

### Steps and overview

* * * * *
To reproduce generated results and figures in the article, nearly all scripts provided in this reposotory.

The scripts provided in the folders from 01 to 10 and data or example input files in the 00_data folder.

#### 1) Create phenoclusters

List of script and their functions for the Step 1:

| **Script** | **Function** |
| --- | --- |
| 00 | [Compare hclust methods](01_clustering/00_compare_hclust_methods.R) |
| 01 | [Create phenoclusters and visualize](01_clustering/01_phenocluster.R) |
| 02 | [Runner script for script 00 and 01](01_clustering/02_LNcluster.bsub) |

</p>

The data for somatic mutation patterns [(SData1.txt)](00_data/SData1.txt) and approved drugs [(SData2.txt)](00_data/SData2.txt)  were accessed on [cBioPortal](https://www.cbioportal.org/) and the [Open Targets Platform](https://platform.opentargets.org/), respectively.

For somatic mutations, each LN phenotype and somatically mutated genes were downloaded from cBioPortal. The LN subtypes-mutated gene binary matrix (0,1) was created based on detected mutated genes without their frequency or position, and the data folder created LN-mutated genes matrix provided as ***SData1.txt***.

For the approved drugs, each phenotype was searched on the Open Targets platform, and each drug (at least phase 3 & 4) was manually searched on public databases to check FDA or EMA approval for clinical usage. The LN subtypes-drug binary matrix (0,1) was created based on this information, and the data folder created LN-drug matrix provided as ***SData2.txt***.

After creating these two data files, hierarchical clustering (hclust) analysis methods were compared using the [**00_Compare_hclust_method.R**](01_clustering/00_compare_hclust_methods.R) script. The script compares hierarchical clustering (hclust) algorithms by generating correlation plots and calculating the Fowlkes-Mallows Index. The resulting dendrograms are analyzed for cophenetic correlation (Pearson correlation coefficient) and Fowlkes-Mallows Index for different cluster counts (k=3, 4, 5). Finally, the script visualizes these correlations using correlation plots saved as TIFF images.

> [!NOTE]
> For somatic mutation data method comparison done with nearly 10 000 genes that are mutated in more than 20% of subtypes becuase of computational restrictions. But, in the next step whole data set used.

After running the **00_Compare_hclust_method.R** script, the selected hclust method is used in the [**01_LNcluster.R**](01_clustering/01_phenocluster.R) script to generate and visualize phenoclusters. The script creates dendrograms using Ward's method and the Jaccard similarity coefficient (in R dist(method = "binary")). The dendrograms are then used to create heatmaps (heatmap.2 from the gplots package) that visualize relationships between drugs or genes and LN subtypes. The resulting heatmaps are saved as TIFF images (drug_plot.tiff and somatic_plot.tiff).

To run these scripts, use the following commands:

```
Rscript --slave --no-restore --no-save 01_clustering/00_compare_hclust_methods.R
Rscript --slave --no-restore --no-save 01_clustering/01_phenocluster.R
```
You can run these commands in your R terminal (not console). If you have access to any HPC, you can submit these scripts with a job runner script. The script 03_LNcluster.bsub is created for IBM LSF job scheduler. This script can easily be converted to commonly used schedulers such as SLURM or PBS. An example SLURM script is provided below:

```
#!/bin/bash
#SBATCH --partition=long                    # Set your partition
#SBATCH --job-name=LNcluster                # Name of the job
#SBATCH --ntasks=1                          # Request a core
#SBATCH --mem=32G                           # Request a total memory limit of 32GB for the job
#SBATCH --output=/your/out/dir/LN_cluster.out  # Give path for the out file
#SBATCH --error=/your/err/dir/LN_cluster.err   # Give path for the err file

# Load modules
module load R/4.3.0

# Set script dir
cd /scripts

# Rscript
Rscript --slave --no-restore --no-save 01_clustering/00_compare_hclust_methods.R
Rscript --slave --no-restore --no-save 01_clustering/01_phenocluster.R

```

* * * * *

#### 2) UKB-QC genetic data and GWAS with regenie

List of script and their functions for the Step 2:

| **Script** | **Function** |
| --- | --- |
| 03 | [QC genotype for regenie step1](02_ukb_gwas/03_qc.sh) |
| 04 | [REGENIE STEP1](02_ukb_gwas/04_regenie_step1.bsub) |
| 05 | [REGENIE STEP2](02_ukb_gwas/05_regenie_step2.bsub) |
| 06 | [Merge REGENIE OUTPUTS](02_ukb_gwas/06_merge_regenie_outputs.sh) |
| 07 | [Filter and format SUMSTATS](02_ukb_gwas/07_filter_format_sumstats.bsub) |

</p>

In this step, provided scripts used to genearate GWAS results from UKB. REGENIE uses a two-step approach. In the first step, original non-imputed genotype data is used, filtering only high-quality genotyped variants: minor allele frequency (MAF) > 1%, minor allele count (MAC) > 5, genotyping rate > 99%, Hardy-Weinberg equilibrium (HWE) test P > 1E−08, <1% missingness. The quality control of genotype data and filtering was done using plink2 software. The script 04_qc.sh combines all chromosomes and creates a list of SNPs that meet QC criteria.

```bash
bash 02_ukb_gwas/03_qc.sh

```

After the QC step, REGENIE step 1 can be run:

```bash
bsub < 02_ukb_gwas/04_regenie_step1.bsub -R "rusage[mem=32G]"

```

REGENIE step 1 will generate a ukb_step1_LM_pred.list file. Using this file and imputed genotype files (bgen), step 2 can be run:

```bash
bsub < 02_ukb_gwas/05_regenie_step2.bsub -R "rusage[mem=32G]"

```

Step 2 will generate GWAS results for each phenotype separated by chromosomes. These files need to be merged by phenotype. To do this, run the merger script:

```bash
bash 02_ukb_gwas/06_merge_regenie_outputs.sh

```

Finally, GWAS summary statistics can be filtered and formatted. The script 08_filter_format_sumstats.bsub creates unique SNP IDs (CHR:POS:REF), converts -log10P to P, and filters Info_score and SPA correction failed tests.

```bash
bsub < 02_ukb_gwas/07_filter_format_sumstats.bsub -R "rusage[mem=32G]"

```

  * * * * *

#### 3) Meta-analysis with METAL


| **Script** | **Function** |
| --- | --- |
| 08 | [METAL parameter](03_meta-analysis/08_s_metalLM.sh) |
| 09 | [Meta-analysis with METAL](03_meta-analysis/10_Metal.bsub) |

</p>

To combine results from UKB, FinnGen and MVP, we performed a fixed-effect meta-analysis using METAL.

```bash
bsub < 03_meta-analysis/10_Metal.bsub -R "rusage[mem=8G]"

```

* * * * *

#### 4) ASSET

List of script and their functions for the Step 3:

| **Script** | **Function** |
| --- | --- |
| 10 | [Parallel ASSET](04_asset/10_asset_parallel.R) |


</p>

To identify pleiotropic variants in a hypothesis-free manner, we conducted an 'Association analysis based on the SubSETs' approach, called ASSET. ASSET is a collection of statistical methods designed to combine association signals from multiple studies or traits, especially when effects are present in only some studies and may be in opposite directions. This tool searches through all potential subsets of studies, adjusts for multiple testing, and identifies the most significant subset contributing to the overall association, while accounting for correlations due to overlapping participants. We ran the ASSET analysis using a custom R script, which enables parallel computation:

```
Rscript --slave --no-restore --no-save 04_asset/10_asset_parallel.R
```


* * * * *


#### 5) LDSC

| **Script** | **Function** |
| --- | --- |
| 11 | [GWASlab](05_ldsc/11_gwaslab_ldsc.bsub) |
| 12 | [LDSC Munge](05_ldsc/12_ldsc_munge.bsub) |
| 13 | [LDSC](05_ldsc//13_ldsc.bsub) |
| 14 | [Combine LDSC Results](05_ldsc//14_combine_ldsc_results.sh) |

</p>

To understand the genetic correlation between LN subtypes and the created phenoclusters, we  performed linkage disequilibrium score regression (LDSC) using LDSC v1.0.1 software. The GWAS summary statistics for LN subtypes and phenoclusters were formatted using the munge_sumstats Python script from LDSC to estimate genetic correlation with HapMap3 variants, as recommended. For the genetic correlation estimation, SNPs with a minor allele frequency (MAF) of less than 5% were excluded, and the MHC region (chr6: 25-35 Mb) was also excluded from this analysis. We downloaded the pre-computed linkage disequilibrium (LD) scores for European ancestry from the Alkes Group website [here](https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays/).

a) GWASlab to format sumstats

```bash
bsub < scripts/11_gwaslab_ldsc.bsub -R "rusage[mem=16G]"

```

b) LDSC munge to format

```bash
bsub < scripts/12_ldsc_munge.bsub -R "rusage[mem=8G]"

```

c) LDSC

```bash
bsub < scripts/13_ldsc.bsub -R "rusage[mem=8G]"

```

d) Combine results

```bash
bash 14_combine_ldsc_results.sh
```

* * * * *


#### 6) HyPrColoc

| **Script** | **Function** |
| --- | --- |
| 11 | [GWASlab](05_ldsc/11_gwaslab_ldsc.bsub) |
| 12 | [LDSC Munge](05_ldsc/12_ldsc_munge.bsub) |
| 13 | [LDSC](05_ldsc//13_ldsc.bsub) |
| 14 | [Combine LDSC Results](05_ldsc//14_combine_ldsc_results.sh) |

</p>


* * * * *

#### 7) susie


* * * * *

#### 8) FLAMES


* * * * *

#### 9) Plots


* * * * *

#### 10) Others

I listed other scripts that we used in various steps to format input/output files.



* * * * *

> [!NOTE]
> This repository is created to enhance reproducibility. It will not be actively maintained.

> [!TIP]
> If you need assistance at any step, please feel free to [email](mailto:murat.guler@dkfz.de) me.

> [!IMPORTANT]
>  It is not possible to share the full datasets used in this pipeline, such as UKBB data. Additionally, some steps require manual verification and the creation of data files.

> [!WARNING]
> This pipeline, or any part of it, should not be used for medical or diagnostic purposes.

> [!CAUTION]
> If you use any part of this pipeline in your work, please make sure to cite the original work software/package and if possible our work.
