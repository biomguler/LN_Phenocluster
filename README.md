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

### Required R Packages

| **R Package** | **Version** | **Link/Citation** |
| --- | --- | --- |
| ASSSET | 2.18.0 | [ASSSET](https://bioconductor.org/packages/release/bioc/html/ASSET.html) |
| tidyverse | 2.0.0 | [tidyverse](https://www.tidyverse.org/) |
| ggplot2 | 3.1.3.1 | [gplots](https://github.com/talgalili/gplots) |
| dendextend | 1.17.1 | [dendextend](https://academic.oup.com/bioinformatics/article/31/22/3718/240978) |
| data.table | 1.14.10 | [data.table](https://github.com/Rdatatable/data.table) |
| parallel | 4.3.0 | [parallel](https://www.R-project.org/) |
| magrittr | 2.0.3 | [magrittr](https://CRAN.R-project.org/package=magrittr) |

### Steps

* * * * *

#### 1) Create phenoclusters

<p style="text-align: justify;"> The data for somatic mutation patterns (SData2.txt) and approved drugs (SData2.txt) were accessed on cBioPortal and the Open Targets Platform, respectively.

For the approved drugs, each phenotype was searched on the Open Targets platform, and each drug (at least phase 3&4) was manually searched on public databases to check FDA or EMA approval for clinical usage. The LN subtypes-drug binary matrix (0,1) was created based on this information, and the data folder created LN-drug matrix provided as ***SData2.txt***.

For the somatic mutations, each LN phenotype and somatically mutated genes were downloaded from cBioPortal. The LN subtypes-mutated gene binary matrix (0,1) was created based on detected mutated genes without their frequency or position, and the data folder created LN-mutated genes matrix provided as ***SData2.txt***.

After creating this two data files, for each dataset, hierarchical clustering (hclust) analysis methods compared by using **00_Compare_hclust_method.R** script. The script aims to compare hierarchical clustering (hclust) algorithms by generating correlation plots and calculating the Fowlkes-Mallows Index. The resulting dendrograms are stored and analyzed for cophenetic correlation (Pearson correlation coefficient) and Fowlkes-Mallows Index for different cluster counts (k=3, 4, 5). Finally, the script visualizes these correlations using correlation plots saved as TIFF images, providing a comprehensive comparison of the clustering methods.

> [!NOTE]
> For somatic mutation data method comparison done with nearly 10 000 genes that are mutated in more than 20% of subtypes becuase of computational restrictions. But, in the next step whole data set used.

After running **00_Compare_hclust_method.R** script, the selected hclust method used in **01_LNcluster.R** script to generate and visualize phenocluster. The script creates dendrograms by using Ward's method and the Jaccard similarity coefficient (in R dist(method = "binary")). The dendrograms are then used to create heatmaps (heatmap.2 from gplots package) that visualize relationships between drugs or genes and LN subtypes. The resulting heatmaps are saved as TIFF images (drug_plot.tiff and somatic_plot.tiff).

</p>

* * * * *

#### 2) QC genetic data and GWAS with regenie

  will be described here full!!

  * * * * *

#### 3) ASSET

  will be described here full!!

#### 4) Meta-analysis with METAL

  will be described here full!!

#### 5) Create FUMA inputs

  will be described here full!!

#### 6) GPS-GEV Test

  will be described here full!!

#### 7) Additional: LDSC

  will be described here full!!

#### 8) 

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
