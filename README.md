LN Phenocluster Manuscript Pipeline
-----------------------------

**Project Name:** Lymphoid neoplasms (LN) phenocluster

**Description:** This repository contains the scripts and pipeline used in the LN phenocluster manuscript. Please follow the scripts listed below to generate analysis results described in the manuscript.

**Citation:** "M. Guler & F. Canzian, Clustering Lymphoid Neoplasms by Somatic Mutation and Drug Usage Profiles: A Multi-trait Genome-wide Association Study, ***submitted***, 2024"

* * * * *

### Required Softwares

| **Software** | **Version** | **Link/Citation** |
| --- | --- | --- |
| regenie | 3.2.1 | [regenie](https://github.com/rgcgithub/regenie) |
| plink | 1.90b6.21 | [plink](https://www.cog-genomics.org/plink/) |
| plink2 | 2.00a3LM | [plink2](https://www.cog-genomics.org/plink/2.0/) |
| gps_cpp | beta | [gps_cpp](https://github.com/twillis209/gps_cpp/tree/remove_po_dependency) |
| R | 4.3.0 or above | [R](https://www.r-project.org/) |

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

### Steps

* * * * *

#### 1) Create phenoclusters

<p style="text-align: justify;"> The data for somatic mutation patterns (SData1.txt) and approved drugs (SData2.txt) were accessed on cBioPortal and the Open Targets Platform, respectively.

For the approved drugs, each phenotype was searched on the Open Targets platform, and each drug (at least phase 3&4) was manually searched on public databases to check FDA or EMA approval for clinical usage. The LN subtypes-drug binary matrix (0,1) was created based on this information, and the data folder created LN-drug matrix provided as ***SData2.txt***.

For the somatic mutations, each LN phenotype and somatically mutated genes were downloaded from cBioPortal. The LN subtypes-mutated gene binary matrix (0,1) was created based on detected mutated genes without their frequency or position, and the data folder created LN-mutated genes matrix provided as ***SData1.txt***.

After creating this two data files, for each dataset, hierarchical clustering (hclust) analysis methods compared by using **00_Compare_hclust_method.R** script. The script aims to compare hierarchical clustering (hclust) algorithms by generating correlation plots and calculating the Fowlkes-Mallows Index. The resulting dendrograms are stored and analyzed for cophenetic correlation (Pearson correlation coefficient) and Fowlkes-Mallows Index for different cluster counts (k=3, 4, 5). Finally, the script visualizes these correlations using correlation plots saved as TIFF images, providing a comprehensive comparison of the clustering methods.

> [!NOTE]
> For somatic mutation data method comparison done with nearly 10 000 genes that are mutated in more than 20% of subtypes becuase of computational restrictions. But, in the next step whole data set used.

After running **00_Compare_hclust_method.R** script, the selected hclust method used in **01_LNcluster.R** script to generate and visualize phenocluster. The script creates dendrograms by using Ward's method and the Jaccard similarity coefficient (in R dist(method = "binary")). The dendrograms are then used to create heatmaps (heatmap.2 from gplots package) that visualize relationships between drugs or genes and LN subtypes. The resulting heatmaps are saved as TIFF images (drug_plot.tiff and somatic_plot.tiff).

Finaly, to run these scripts please run code below:

```r
Rscript --slave --no-restore --no-save scripts/00_compare_hclust_methods.R
Rscript --slave --no-restore --no-save scripts/01_phenocluster.R
```
You can run these two one-line codes in you R terminal (**not console**), and if you have access to any HPC, you can push these script with a job runner script. The script **03_LNcluster.bsub** is created for IBM LSF job scheduler.
But, this script easly can be converted to commonly used scheduler such as SLURM or PBS. In pipeline all runner scripts were created for LSF and as an example how to convert this scripts to SLURM or any other platform an example code snippet for SLURM given in below:

```bash
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
Rscript --slave --no-restore --no-save scripts/00_compare_hclust_methods.R
Rscript --slave --no-restore --no-save scripts/01_phenocluster.R

```

List of script and their functions for the step1:

| **Script** | **Function** |
| --- | --- |
| 00 | [Compare hclust methods](scripts/00_compare_hclust_methods.R) |
| 01 | [Create phenoclusters and visualize](scripts/01_phenocluster.R) |
| 03 | [Runner script for script 00 and 01](scripts/03_LNcluster.bsub) |

</p>

* * * * *

#### 2) QC genetic data and GWAS with regenie

REGENIE uses a two-step approach. In the first step, original non-imputed genotype data is used, filtering only high-quality genotyped variants: minor allele frequency (MAF) > 1%, minor allele count (MAC) > 5, genotyping rate >99%, Hardy-Weinberg equilibrium (HWE) test P > 1E−08, <1% missingness. The quality control of genotype data and filtering was done by using plink2 software. The script **04_qc.sh** combine all chromosemes and create list of SNPs which meet QC criteria. 

```bash
bash scripts/04_qc.sh

```

After the QC step, we can run REGENIE step 1.

```bash
bsub < scripts/05_regenie_step1.bsub -R "rusage[mem=32G]"

```

The REGENIE step 1 will generate **ukb_step1_LM_pred.list** file and by using this and imputated genotype files (bgen), we can run step 2.

```bash
bsub < scripts/06_regenie_step2.bsub -R "rusage[mem=32G]"

```

The step 2, will generate GWAS results for each phenotype seperated by chromosomes. We need to merge these files by phenotype. To do this we can run merger script.

```bash
bash scripts/07_merge_regenie_outputs.sh

```

Finally, we can filter and format the GWAS summary statistics. The script **08_filter_format_sumstats.bsub**, creates unique SNP ids (CHR:POS:REF:TEST), covert -log10P to P and filter Info_score and SPA correction failed tests. 

```bash
bsub < scripts/08_filter_format_sumstats.bsub -R "rusage[mem=32G]"

```

List of script and their functions for the step2:

| **Script** | **Function** |
| --- | --- |
| 04 | [QC genotype for regenie step1](scripts/04_qc.sh) |
| 05 | [REGENIE STEP1](scripts/05_regenie_step1.bsub) |
| 06 | [REGENIE STEP2](scripts/06_regenie_step2.bsub) |
| 07 | [Merge REGENIE OUTPUTS](scripts/07_merge_regenie_outputs.sh) |
| 08 | [Filter and format SUMSTATS](scripts/08_filter_format_sumstats.bsub) |

</p>


  * * * * *

#### 3) ASSET and Meta-analysis with METAL

As a hypothesis-free approach to identify pleiotropic variants, we conduct Association analysis based on SubSETs approach called ASSET [7]. ASSET is a collection of statistical methods tailored to combine association signals from multiple studies or traits, particularly when effects are present in only some studies and may be in opposite directions. The tool searches through all potential subsets of studies, adjusts for multiple testing, and identifies the most significant subset contributing to the overall association, accounting for correlations due to overlapping participants. We ran ASSET analysis with a custom R script which enables to compute the test parallel computation. 

```r
Rscript --slave --no-restore --no-save scripts/09_asset_parallel.R
```
To compare ASSET, phenocluster and traditional meta-analysis we ran meta-analysis with METAL. Because we have 7 phenocluster, we create seven script for METAL, 10_s_metalLM(**1-7**).sh, and we run all this script with **10_Metal.bsub**.

```bash
bsub < scripts/10_Metal.bsub -R "rusage[mem=8G]"

```

List of script and their functions for the step3:

| **Script** | **Function** |
| --- | --- |
| 09 | [Parallel ASSET](scripts/09_asset_parallel.R) |
| 10 | [Meta-analysis with METAL](scripts/10_Metal.bsub) |

</p>



#### 4) GPS-GEV Test

  will be described here full!!

#### 5) Create FUMA inputs

  will be described here full!!

#### 6) Additional: LDSC

  will be described here full!!

#### 7) Plots

  will be described here full!!


* * * * *

> [!NOTE]
> This repository is created to enhance reproduciblity. It will not actively maintained.

> [!TIP]
> If you need any help at any step please send me an [email](mailto:murat.guler@dkfz.de).

> [!IMPORTANT]
> It is impossible to share full data sets that are used in the pipeline like UKBB data and some steps need manual check and create data files.

> [!WARNING]
> The pipeline or any part of work can not use medical/diagnostic proposes.

> [!CAUTION]
> If you use any part of this pipeline in your work please cite our orginal work!
