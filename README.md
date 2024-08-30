LN Phenocluster Manuscript Pipeline
-----------------------------

**Project Name:** Lymphoid neoplasms (LN) phenocluster

**Description:** This repository contains the scripts and pipeline used in the LN phenocluster manuscript. Please follow the scripts listed below to generate the analysis results described in the manuscript.

**Citation:** "M. Guler & F. Canzian, Clustering of Lymphoid Neoplasms by Cell of Origin, Somatic Mutation and Drug Usage Profiles: A Multi-trait Genome-Wide Association Study, ***submitted***, 2024"

* * * * *

### Required Softwares

| **Software** | **Version** | **Link/Citation** |
| --- | --- | --- |
| regenie | 3.2.1 | [regenie](https://github.com/rgcgithub/regenie) |
| plink | 1.90b6.21 | [plink](https://www.cog-genomics.org/plink/) |
| plink2 | 2.00a3LM | [plink2](https://www.cog-genomics.org/plink/2.0/) |
| gps_cpp | beta | [gps_cpp](https://github.com/twillis209/gps_cpp/tree/remove_po_dependency) |
| R | 4.3.0 or above | [R](https://www.r-project.org/) |
| METAL | 2020-05-05 | [METAL](https://github.com/statgen/METAL) |
| LDSC | v1.0.1 | [LDSC](https://github.com/bulik/ldsc) |
| GWASLab | v3.4.31 | [GWASLab](https://github.com/Cloufield/gwaslab) |

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

The data for somatic mutation patterns (SData1.txt) and approved drugs (SData2.txt) were accessed on cBioPortal and the Open Targets Platform, respectively.

For the approved drugs, each phenotype was searched on the Open Targets platform, and each drug (at least phase 3 & 4) was manually searched on public databases to check FDA or EMA approval for clinical usage. The LN subtypes-drug binary matrix (0,1) was created based on this information, and the data folder created LN-drug matrix provided as ***SData2.txt***.

For somatic mutations, each LN phenotype and somatically mutated genes were downloaded from cBioPortal. The LN subtypes-mutated gene binary matrix (0,1) was created based on detected mutated genes without their frequency or position, and the data folder created LN-mutated genes matrix provided as ***SData1.txt***.

After creating these two data files, hierarchical clustering (hclust) analysis methods were compared using the **00_Compare_hclust_method.R** script. The script compares hierarchical clustering (hclust) algorithms by generating correlation plots and calculating the Fowlkes-Mallows Index. The resulting dendrograms are analyzed for cophenetic correlation (Pearson correlation coefficient) and Fowlkes-Mallows Index for different cluster counts (k=3, 4, 5). Finally, the script visualizes these correlations using correlation plots saved as TIFF images.

> [!NOTE]
> For somatic mutation data method comparison done with nearly 10 000 genes that are mutated in more than 20% of subtypes becuase of computational restrictions. But, in the next step whole data set used.

After running the **00_Compare_hclust_method.R** script, the selected hclust method is used in the **01_LNcluster.R** script to generate and visualize phenoclusters. The script creates dendrograms using Ward's method and the Jaccard similarity coefficient (in R dist(method = "binary")). The dendrograms are then used to create heatmaps (heatmap.2 from the gplots package) that visualize relationships between drugs or genes and LN subtypes. The resulting heatmaps are saved as TIFF images (drug_plot.tiff and somatic_plot.tiff).

To run these scripts, use the following commands:

```
Rscript --slave --no-restore --no-save scripts/00_compare_hclust_methods.R
Rscript --slave --no-restore --no-save scripts/01_phenocluster.R
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
Rscript --slave --no-restore --no-save scripts/00_compare_hclust_methods.R
Rscript --slave --no-restore --no-save scripts/01_phenocluster.R

```

List of script and their functions for the Step 1:

| **Script** | **Function** |
| --- | --- |
| 00 | [Compare hclust methods](scripts/00_compare_hclust_methods.R) |
| 01 | [Create phenoclusters and visualize](scripts/01_phenocluster.R) |
| 03 | [Runner script for script 00 and 01](scripts/03_LNcluster.bsub) |

</p>

* * * * *

#### 2) QC genetic data and GWAS with regenie

REGENIE uses a two-step approach. In the first step, original non-imputed genotype data is used, filtering only high-quality genotyped variants: minor allele frequency (MAF) > 1%, minor allele count (MAC) > 5, genotyping rate > 99%, Hardy-Weinberg equilibrium (HWE) test P > 1E−08, <1% missingness. The quality control of genotype data and filtering was done using plink2 software. The script 04_qc.sh combines all chromosomes and creates a list of SNPs that meet QC criteria.

```bash
bash scripts/04_qc.sh

```

After the QC step, REGENIE step 1 can be run:

```bash
bsub < scripts/05_regenie_step1.bsub -R "rusage[mem=32G]"

```

REGENIE step 1 will generate a ukb_step1_LM_pred.list file. Using this file and imputed genotype files (bgen), step 2 can be run:

```bash
bsub < scripts/06_regenie_step2.bsub -R "rusage[mem=32G]"

```

Step 2 will generate GWAS results for each phenotype separated by chromosomes. These files need to be merged by phenotype. To do this, run the merger script:

```bash
bash scripts/07_merge_regenie_outputs.sh

```

Finally, GWAS summary statistics can be filtered and formatted. The script 08_filter_format_sumstats.bsub creates unique SNP IDs (CHR:POS:REF), converts -log10P to P, and filters Info_score and SPA correction failed tests.

```bash
bsub < scripts/08_filter_format_sumstats.bsub -R "rusage[mem=32G]"

```

List of script and their functions for the Step 2:

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

To identify pleiotropic variants in a hypothesis-free manner, we conducted an Association analysis based on the SubSETs approach, called ASSET. ASSET is a collection of statistical methods designed to combine association signals from multiple studies or traits, especially when effects are present in only some studies and may be in opposite directions. This tool searches through all potential subsets of studies, adjusts for multiple testing, and identifies the most significant subset contributing to the overall association, while accounting for correlations due to overlapping participants. We ran the ASSET analysis using a custom R script, which enables parallel computation:

```
Rscript --slave --no-restore --no-save scripts/09_asset_parallel.R
```
To compare the results from ASSET, phenocluster, and traditional meta-analysis, we performed a meta-analysis using METAL. Since we have seven phenoclusters, we created seven separate scripts for METAL (10_s_metalLM(1-7).sh), and executed all these scripts with 10_Metal.bsub:

```bash
bsub < scripts/10_Metal.bsub -R "rusage[mem=8G]"

```

List of script and their functions for the Step 3:

| **Script** | **Function** |
| --- | --- |
| 09 | [Parallel ASSET](scripts/09_asset_parallel.R) |
| 10 | [Meta-analysis with METAL](scripts/10_Metal.bsub) |

#### 4) GPS-GEV Test and LDSC

To understand genetic correlation between LN subtypes and created phenocluster, we first performed linkage disequilibrium score regression (LDSC) by using LDSC v1.0.1 software. The GWAS summary statistics of LN subtypes and phenoclusters were formatted with munge_sumstats python script from LDSC to estimate genetic correlation with HapMap3 variants as recommended. For genetic correlation estimation, SNPs were excluded if the MAF was smaller than 5% and the MHC region (chr6: 25-35 Mb) was excluded from this analysis. The pre-computed linkage disequilibrium (LD) scores for European ancestry were downloaded from the Alkes Group website (https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays/).

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

Because our effective sample size [Neff=4 /(1/Ncases+1/Ncontrols)] was relatively small, the majority of pairwise correlation tests failed. Therefore, we used an alternative method called the genome-wide pairwise-association signal sharing (GPS) test [9]. As a non-parametric test, the test does not provide quantitative measurement of correlation like LDSC, but provides pairwise independence. Instead it is only evidence against a null hypothesis of bivariate independence. 
Recently, the GPS was enhanced by fitting it with a generalized extreme value distribution instead of using the standard exponential transformation as initially proposed. We used GPS-GEV test. We only removed the MHC region (chr6 29624758-33170276). We used computeGpsCLI application to compute the GPS test statistic for P values from a pair of GWAS and permuteTraitsCLI application used to generate null realizations of the GPS test statistic with 100 permutation. Finally, P values GPS statistic obtained with fit_gevd_and_compute_pvalue.R script which fits a generalized extreme value distribution (GEVD) to null realizations of the GPS statistic and reports a P value.

a) To compute GPS:
```bash
bsub < scripts/15_computeGPScli.bsub -R "rusage[mem=200G]"

```

b) To permute GPS:
```bash
bsub < scripts/16_permuteTraitsCLI.bsub -R "rusage[mem=200G]"

```

c) To fit GEVD and claculate P value:
Note: Script 18_compute_pvalue.bsub uses **17_fit_gevd_compute_pvalue.R**.
```bash
bsub < scripts/18_compute_pvalue.bsub -R "rusage[mem=200G]"

```

| **Script** | **Function** |
| --- | --- |
| LDSC |  |
| 11 | [GWASlab](scripts/11_gwaslab_ldsc.bsub) |
| 12 | [LDSC Munge](scripts/12_ldsc_munge.bsub) |
| 13 | [LDSC](scripts/13_ldsc.bsub) |
| 14 | [Combine LDSC Results](scripts/14_combine_ldsc_results.sh) |
| GPS-GEV |  |
| 15 | [Compute GPS](scripts/15_computeGPScli.bsub) |
| 16 | [Permute GPS](scripts/16_permuteTraitsCLI.bsub) |
| 17 | [Fit GEV-Pvalue R script ](scripts/17_fit_gevd_compute_pvalue.R) |
| 18 | [Runner for script 17](scripts/18_compute_pvalue.bsub) |




#### 5) Create FUMA inputs

FUMA requests certain GWAs summary statistics format and less than 600Mb size. To achive this, we used different custom R scripts for Regenie, METAL and ASSET GWAS summary statistics. 

After FUMA results for each GWAS summary statistics, results need to be combined and summurized the script **22_fuma2functional.R** used for this purpose.

a) Regenie to FUMA

```
Rscript --slave --no-restore --no-save scripts/19_regenie2fuma_input.R
```

b) METAL to FUMA

```
Rscript --slave --no-restore --no-save scripts/20_metal2fuma.R
```

c) ASSET to FUMA

```
Rscript --slave --no-restore --no-save scripts/21_asset2fuma.R
```

d) FUMA to summary functional tables

```
Rscript --slave --no-restore --no-save scripts/22_fuma2functional.R
```



#### 7) Plots
The result of GPS-GEV and LDSC plotted rogether with using custom R script.

a) To generate corrplot for GPS-GEV and LDSC

```
Rscript --slave --no-restore --no-save scripts/23_gps_corrplot.R
```

b) To generate custom forest plots

Note: The script is adsapted from Katherine Hoffman's blogpage. The [source](https://www.khstats.com/blog/forest-plots/#just-the-code) (the last accessed date: 30/08/2024).

```
Rscript --slave --no-restore --no-save scripts/24_forest_plot.R
```

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
