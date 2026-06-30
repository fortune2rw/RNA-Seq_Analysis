# Bulk RNA-Seq 
##  TCGA RNA-seq Exploratory Analysis (Tumour vs Normal)

## Overview
Exploratory data analysis (EDA) of RNA-seq gene expression data obtained from The Cancer Genome Atlas (TCGA) using the GDC API.

- The aim of the analysis is to perform quality control checks.
- Identify patterns prior to differential expression analysis.
- measure the similarity of samples in terms of quantified expression level across profiles.


## Method
Data was retrieved using the GDC API via the `TCGAbiolinks` package.
### Query parameters:
- **Project**: (TCGA-LUAD)
- **Data type**: Gene Expression Quantification
- **Workflow**: STAR - Counts
- **Sample types**:
  - Primary Tumor (CASE)
  - Solid Tissue Normal (CTRL)
 
---
## Data Cleaning & Transformation
- Ensembl gene IDs were mapped to gene symbols using:
  - `org.Hs.eg.db`
- Duplicate gene were removed using "janitor" package
- Safety checked alighment of gene symbols and values to ensure data integrity

 **Limitation**:  
By removing duplicates, biological releveant information may be lost since multiple Ensembl IDs can represent different transcript isoforms of the same gene.

**Future Works**:
Researchers may aggregate this duplicates( e.g.  mean or sum) rather than totally get rid of them

## Normalization:
Before comparing gene expression across samples, several technical biases were considered:
### sources of bias:
1. GC content bias
2. Longer genes will have longer number of reads
3. Library composition can be different in two type of biological conditions
4. Library size varies between samples coming from diffeent lanes of flow cell in the sequencing machine.
5. Read coverage of a transcript can be biased and non uniformly distributed along the transcript.

### Normalization approches
1. To address the library size bias(sequencing depth), read counts per gene should be divided by each gene count by a certain value and multiplied by 10^6. These values are referedd to as CPM ( count per million).

- Other metrics that imoroves CPM are RPKM/FPKM(read/fragment per kilo base of million read) and TPM (transcript per million).

## Data Visualisation
### Selected top 50 significant genes ###
**rationale* : Samples can be represntaed by a couple of principal varibales instead of thousands of genes. This is useful for visualization, clustering and predictive modeling.

### Principal Component Analsysi (PCA)
- PCA analysis highlited  any initital quality control issues that may arise from the sample which may affect potential analsysis results.
- PCA also highlight any batch effects across cohorts

### Clustering Heatmap
- Clustering samples by row (gene) using "pheatmap" to find out whcih genes have the same expression level across Control & Cases.
- Column clustering (samples): Groups samples based on transcriptional similarity

### Correlation Heatmap
- Compute correlation score between each pair of samples.
- Visualized sample correlation to identify which samples
**rationale* : Correlation analysis improves clustering quality.


## Differential Expression Analysis:
- DESeq2.
  
## Packages  
``` r
install.packages(c(
  "BiocManager", "ggfortify", "pheatmap", "reshape", "ggplot2", "corrplot", "dplyr",
"stats", "janitor"
))

