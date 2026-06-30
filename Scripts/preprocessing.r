# packages
install.packages("BiocManager")
install.packages("janitor")
install.packages("ggfortify")
install.packages("corrplot")
install.packages("ggrepel")
BiocManager::install(c("TCGAbiolinks","AnnotationDbi", "org.Hs.eg.db",
                       "recount3","DESeq2", "clusterProfiler"), quietly = TRUE)

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
library(DESeq2)
library(reshape)
library(ggrepel)
library(clusterProfiler)

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
LUAD_data_mat <- assay(LUAD_data) #obtain count 
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

#change rowname from original LUAD matrix to new esemblid
rownames(LUAD_data_mat) <- ensembl_ids

# most Esembl id has the same gene symbol,lets verify
gene_symbol %>%
  janitor::get_dupes(SYMBOL)

#remove NA
gene_symbol %>%
  janitor::get_dupes(SYMBOL) %>%
  filter(!is.na(SYMBOL))

  #remove duplicate
gene_symbol_uniq <- gene_symbol %>%
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
  labs(title = "Heatmap showing sample correlation",
       x = NULL, y = NULL, fill = "Correlation")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 10))

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

#create metadata
metadata <- data.frame(
  sample_name = colnames(subset_LUAD_data_mat),
  condition = c(rep("tumour", 10), rep("normal", 10))
)
#convert condition variable to factor
metadata$condition <- as.factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref = "normal")

dds <- DESeqDataSetFromMatrix(countData = subset_LUAD_data_mat,
                              colData = metadata,
                              design = ~ condition
                              )
#For each gene, we count the total number of reads for that gene in all samples 
#and remove those that don't have at least 1 read. 
dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]
#Run DESEQ 
dds <- DESeq(dds)
results <- results(dds, alpha = 0.05)
results_df <- as.data.frame(results)
#remove rows with missing values & get significant genes with < 0.05 p value
sig_gene <- results[!is.na(results$padj) & results$padj < 0.05, ]
#get up-regulated genes
sig_gene <- as.data.frame(sig_gene) ##change to dataframe
upregulated <- sig_gene %>%
  filter(log2FoldChange > 2 )
#get down regulated genes
downregulated <- sig_gene %>%
  filter(log2FoldChange < -2 )

sig_gene_all <- rbind(upregulated, downregulated)
#print output to csv file
write.csv(upregulated, file = "upregulated_gene.csv", row.names = T)
write.csv(downregulated, file = "downregulated_gene.csv", row.names = T)
#diagnostic plots
MA_plot <- plotMA(results, ylim = c(-5, 5))

##Data visualization
volcano_data <- as.data.frame(results) %>% #prepare the data
  filter(!is.na(padj)) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 2  ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -2 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )
volcano_data$gene <- rownames(volcano_data)

# Get top all significant genes  & arrange by smallest p value to label
top_genes <- volcano_data %>%
  filter(significance != "Not Significant") %>%
  arrange(padj)

#volcano plot
ggplot(volcano_data, aes(x = log2FoldChange, 
                         y = -log10(padj), 
                         color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "gray"
  )) +
  geom_text_repel(data = top_genes, 
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 10) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Lung Cancer vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Expression"
  ) +
  theme_minimal()


##Plot Heatmap
# Get top 50 significant genes by padj
top50_genes <- sig_gene %>%
  arrange(padj) %>%
  head(50)
# Extract normalized counts for those genes
normalized_counts <- counts(dds, normalized = TRUE)
heatmap_data <- normalized_counts[rownames(top50_genes), ]

# Log transform for better visualization
heatmap_data <- log2(heatmap_data + 1)

# Create annotation for samples
annotation_col <- data.frame(
  condition = metadata$condition,
  row.names = metadata$sample_name
)

# Plot 
pheatmap(heatmap_data,
         annotation_col = annotation_col,
         annotation_colors = list(
           condition = c(tumour = "red", normal = "blue")
         ),
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap: Top 50 DEGs in Lung Cancer vs Normal")

