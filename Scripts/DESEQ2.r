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

