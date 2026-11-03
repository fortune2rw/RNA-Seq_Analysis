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

###############################################################################
# Differential Expression Analysis 
###############################################################################

#create metadata
# Create metadata pulling groups organically so they scale perfectly
metadata <- data.frame(
  sample_name = colnames(subset_LUAD_data_mat),
  condition   = ifelse(met$definition == "Solid Tissue Normal", "normal", "tumour")
)
write.csv(metadata, "result/DEG/DESEQ_metadata.csv")
print(paste("Saved metadata table to a result subfolder"))

#convert condition variable to factor
metadata$condition <- as.factor(metadata$condition) #converting to factor tells r that it is a categorical group rather than a text character.
metadata$condition <- relevel(metadata$condition, ref = "normal") #when comparing two condition/groups u always need a baseline to compare with. in this scenario normal samples compared to tumour sample.
#we want to know the genes that are diferentially expressed in tumour sample
#we relevel so to tell r which baseline to choose first, be default r picks alphabetically.
# if we dont relevel, the final result is counted backwards perhaps r picks the wrong baseline. upregulated genes in tumour will be intepreted as downregulated

#create a Desq dataset
dds <- DESeqDataSetFromMatrix(countData = subset_LUAD_data_mat,
                              colData = metadata,
                              design = ~ condition
)

#For each gene, we count the total number of reads for that gene in all samples 
#and remove those that don't have at least 1 read. 
dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]

#Run DESEQ 
dds <- DESeq(dds, test = "Wald", sfType = "poscount")
results <- results(dds, alpha = 0.05, pAdjustMethod ="fdr")
results_df <- as.data.frame(results)
write.csv(results_df, file = "result/DEG/DESEQ_result.csv")
print(paste("Saved Deseq table to a result subfolder"))

# create histogram plot of p-values
png("padj_histogram.png", width = 800, height = 600, res = 120)

hist(results$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = " Frequencies of padj-values")
dev.off()

#plot dispersion estimate 
png(filename = "dispersion_estimate.png", width = 1300, height = 840, res = 250)

plotDispEsts(dds, main = "Dispersion Estimates")

dev.off()

#remove rows with missing values & get significant genes with < 0.05 p value
sig_gene <- results[!is.na(results$padj) & results$padj < 0.05, ]

#get up-regulated genes
sig_gene <- as.data.frame(sig_gene) ##change to dataframe
upregulated <- sig_gene %>%
  filter(log2FoldChange > 2 )
write.csv(sig_gene, file = "result/DEG/sig_gene.csv")


#get down regulated genes
downregulated <- sig_gene %>%
  filter(log2FoldChange < -2 )

sig_gene_all <- rbind(upregulated, downregulated)
#print output to csv file
write.csv(upregulated, file = "result/DEG/upregulated_gene.csv", row.names = T)
write.csv(downregulated, file = "result/DEG/downregulated_gene.csv", row.names = T)
#diagnostic plots
png(filename = "MA_plot.png", width = 1300, height = 800, res = 250)
MA_plot <- plotMA(results, ylim = c(-5, 5))
dev.off()

################    DATA VISUALIZATION   ################ ################ 
################ ################  ################ ################ 


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
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, 
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
ggsave("volcanoplot.png", plot = volcano_plot, width = 10, height =12, dpi = 300)


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
         main = "Heatmap: Top 50 DEGs in Tumour vs Normal Sample")


#   General expression data visualization
dat <- log10(counts(dds, normalized = T) + 1) # extract normalized counts

# box-and-whisker plot
#define group string (metadata$condition) as gs
gs <- as.factor(metadata$condition)
lbl <- "log10(raw counts + 1)"
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
boxplot(dat[,ord], boxwex=0.6, notch=T, main="Boxplot", ylab="lg(norm.counts)", outline=F, las=2, col=gs[ord])
legend("topleft", legend = levels(gs), fill = palette()[1:nlevels(gs)], bty = "n")

# UMAP plot (multi-dimensional scaling)
dat <- dat[!duplicated(dat), ] # first remove duplicates
par(mar=c(3,3,2,6), xpd=TRUE, cex.main=1.5)
ump <- umap(t(dat), n_neighbors = 5, random_state = 123)
plot(ump$layout, main="UMAP plot", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset = c(-0.15, 0), legend = levels(gs), pch = 20,
       col = palette()[1:nlevels(gs)], title = "Group", pt.cex = 1.5)

