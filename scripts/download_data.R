# libraries
library(dplyr)
library(SummarizedExperiment)
library(TCGAbiolinks)

# ==============================================================================
# DATA IMPORTATION
# ==============================================================================
#Query to download the data from GDC
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  experimental.strategy = "RNA-Seq",
                  access = "open",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

#download a subset of normal tissue and tumour from query, 10 each
set.seed(20)
samples <- getResults(query)
CTRL <- unique(samples$cases[samples$sample_type == "Solid Tissue Normal"])[1:5]

CASE  <- unique(samples$cases[samples$sample_type == "Primary Tumor" ])[1:5]

selected_samples <- c(CASE, CTRL)

selected_sample_info <- samples[samples$cases %in% selected_samples, ]

# Create a subset query using the barcodes
# NOTE: Adjusted column name to match TCGAbiolinks standard "cases" barcode column
# Study population are selected based on individual(cases) and not by thier sample
query_subset <- GDCquery(project = "TCGA-LUAD",
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts",
                         experimental.strategy = "RNA-Seq",
                         access = "open",
                         barcode = selected_sample_info$sample.submitter_id)

#Download the data
GDCdownload(query_subset, directory = "data/raw")
LUAD_data <- GDCprepare(query_subset, directory = "data/raw") 
saveRDS(LUAD_data, file = "data/raw/.rds")

#Pull out information from data
LUAD_data_mat <- assay(LUAD_data) #obtain count 
met <- colData(LUAD_data)

print("Downloaded TCGA-LUAD cohort, group of 5 CTRL and 5 CASE")
print(paste("Saved matrix count to raw data folder, row:",length(LUAD_data_mat),", column:", length(colnames(LUAD_data_mat))))


