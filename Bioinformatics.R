# packages
install.packages("BiocManager")
install.packages("janitor")
install.packages("DESeq2")
install.packages("ggfortify")
install.packages("corrplot")
BiocManager::install(c("TCGAbiolinks","AnnotationDbi", "org.Hs.eg.db",
                       "recount3"), force = TRUE)

# libraries
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(recount3)
library(readr)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(janitor)
library(pheatmap)
library(stats)
library(ggplot2)
library(ggfortify)
library(corrplot)
library(reshape)

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

samples <- getResults(query)
CTRL <- samples[samples$sample_type == "Solid Tissue Normal",][1:10,]
CASE <- samples[samples$sample_type == "Primary Tumor",][1:10,]
selected_samples <- rbind(CASE, CTRL)

#subset query using barcode
query_subset <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  experimental.strategy = "RNA-Seq",
                  access = "open",
                  barcode = selected_samples$sample.submitter_id
)

#Download the data
GDCdownload(query_subset)
LUAD_data <- GDCprepare(query_subset) 


#merge files
files <- list.files(path = "~/Desktop/Bioinformatics/GDCdata/TCGA-LUAD/Transcriptome_Profiling/Gene_Expression_Quantification", full.names = T)

#Pull out information from data
LUAD_data_mat <- assay(LUAD_data) #count 
met <- colData(LUAD_data)

# ==============================================================================
# DATA CLEANING
# ==============================================================================
ensemblids <-rownames(LUAD_data_mat)

#remove version number 
ensembl_ids <- ensemblids %>%
  stringr::str_replace("\\.[0-9]+$", "")

#mapping gene id to gene symbol
gene_symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = ensembl_ids,
                                     columns = "SYMBOL",
                                     keytype = "ENSEMBL",
                                     multiVar = "first")
#chnage rownmaes from original LUAD matrix to new esemblid
rownames(LUAD_data_mat) <- ensembl_ids

# most Esembl id has the same gene symbol,lets verify
gene_symbol %>%
  janitor::get_dupes(SYMBOL)

#remove NA
gene_symbol %>%
  janitor::get_dupes(SYMBOL) %>%
  filter(!is.na(SYMBOL))

  #remove duplicate
gene_symbol_uniq <-gene_symbol %>%
  filter(!is.na(SYMBOL))%>%
    distinct(SYMBOL, .keep_all = TRUE)%>%
    distinct(ENSEMBL, .keep_all = TRUE)

##subset the original matrix
subset_LUAD_data_mat <- LUAD_data_mat[gene_symbol_uniq$ENSEMBL, ]

all.equal(rownames(subset_LUAD_data_mat),
          gene_symbol_uniq$ENSEMBL) #verify alignment

rownames(subset_LUAD_data_mat) <- gene_symbol_uniq$SYMBOL
head(subset_LUAD_data_mat)

### chnage columns to approrpaite sample 
colnames(subset_LUAD_data_mat) <- met$definition
colnames(subset_LUAD_data_mat)
#Rename columns 
colnames(subset_LUAD_data_mat)[colnames(subset_LUAD_data_mat) == "Solid Tissue Normal"] <- "Control"
colnames(subset_LUAD_data_mat)[colnames(subset_LUAD_data_mat) == "Primary solid Tumor"] <- "Case"
colnames(subset_LUAD_data_mat) <- make.unique(colnames(subset_LUAD_data_mat))


###############################################################################
# Exploratory Data Analysis
###############################################################################

##Select top 50 significant genes 
# compute the variance of each gene across samples
genes_V <- apply(subset_LUAD_data_mat, 1, var)
# sort by decreasing order and select top 100
selected_genes <- names(sort(genes_V, decreasing = TRUE))[1:50]

#subset selected genes from original matrix
gene_mat <- subset_LUAD_data_mat[selected_genes,]

#change column names 
colnames(gene_mat) <- c(paste0("Case_", 1:10),
paste0("Control_", 1:10))

##plot heatmap
annotation_col <- data.frame(
  Group = ifelse(grepl("^Control", colnames(gene_mat)), "CTRL", "CASE")
)
rownames(annotation_col) <- colnames(gene_mat) #clustering samples with anaotation

pheatmap(gene_mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = TRUE,
         annotation_col = annotation_col,
         main = "Heatmap Clustered by genes", fontsize_number = 15)


##PCA to see clustering of replicates as a scatter plot in 2 dimension
    
#transpose the matrix
M <- t(subset_LUAD_data_mat[selected_genes,])

#transform counts to log2scale
M <- log2(M + 1)

#compute PCA 
pcaResult <- prcomp(M, scale. = T)
pca_df <- data.frame(
  PC1 = pcaResult$x[,1],
  PC2 = pcaResult$x[,2],
  Group = ifelse(grepl("^Control", rownames(M)), "CTRL", "CASE")
)

#plot
autoplot(pcaResult, data = pca_df, color = "Group", size = 4)+
  labs( title = "PCA Scatter Plot  in 2D")
summary(pcaResult)

## Correlation plot
cor_matrix <- cor(gene_mat)
cor_long <- melt(cor_matrix)

ggplot(cor_long,
       aes(x = X1 , y = X2, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
  labs( title = "Correlation PLot shoing correlation coefficeinct as percentage (%)")+
  theme_minimal()

#plot
corrplot(cor_matrix, addCoefasPercent = TRUE, addCoef.col = "white", tl.srt = 25)

#pairwise correlation sample displayed as heatmap
pheatmap(mat = cor_matrix,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         annotation_col = annotation_col,
         main = "Pairwise Correlation Sample Displayed as Heatmap")



###############################################################################
# Differential Expression Analysis 
###############################################################################

