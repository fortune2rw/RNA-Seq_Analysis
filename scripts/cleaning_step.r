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
library(umap)

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

### chnage column names mapped by metadata
clean_groups <- ifelse(met$definition == "Solid Tissue Normal", "Control", "Case")
# Assign these clean, unique names directly to your matrix
colnames(subset_LUAD_data_mat) <- make.unique(clean_groups) 
print(colnames(subset_LUAD_data_mat)) #the matrix had 11 columns instead of 10. One of the patient( must have had more than one sample downloaded)


###############################################################################
# Exploratory Data Analysis
###############################################################################

##Select top 50 significant genes 
# compute the variance of each gene across samples
# --- Top 50 most variable genes ---
# log-transform first so highly expressed genes don't dominate the variance
log_mat <- log2(subset_LUAD_data_mat + 1)

genes_V <- apply(log_mat, 1, var)
# sort by decreasing variance and select top 50
selected_genes <- names(sort(genes_V, decreasing = TRUE))[1:50]

# subset selected genes (log2 scale)
gene_mat <- log_mat[selected_genes, ]

# --- Rename columns from the sample metadata, not column position ---
# met$definition is in the same order as the matrix columns
group <- ifelse(met$definition == "Solid Tissue Normal", "Control", "Case")
stopifnot(length(group) == ncol(gene_mat))
table(group)   # sanity check: should show 5 Case, 5 Control

colnames(gene_mat) <- paste0(group, "_", ave(seq_along(group), group, FUN = seq_along))


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


