###############################################################################
# GO Enrichment Analysis 
###############################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

# in this analysis, we will match gene to their biological process, molecular process and cellular process
#future work: pathway enrichment
#read deseq result from path 
clust_df <- read.csv("result/DEG/DESEQ_result.csv")

#define sign genes 
sig_gene_clust <- clust_df$X[!is.na(clust_df$padj) & clust_df$padj < 0.05]
background_genes <- clust_df$X

# first we need to convert gene symbol to entrez id so cluster profiler can read it
# we can use the bitr function 
gene_set <-clusterProfiler::bitr(geneID = sig_gene_clust,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db"
                      )

#turn background gene symbol into ENTREZID , same as above step for sig gene
backgrund_set <-clusterProfiler::bitr(geneID = background_genes,
                                 fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = "org.Hs.eg.db"
)

#compute gene set enrichments
enrich_GO <- clusterProfiler::enrichGO(gene = gene_set$ENTREZID,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pvalueCutoff = 0.01,
                          pAdjustMethod = "fdr",
                          OrgDb = "org.Hs.eg.db",
                          universe = backgrund_set$ENTREZID
                          )
enrich_GO_df <- as.data.frame(enrich_GO)
write.csv(enrich_GO_df, file = "result/Geneset_enrichment/enrich_GO_BiologicalProcess_results.csv", row.names = FALSE)
print("Successfully saved biological process gene enrichment to result folder")

#Quick visalisation
png("biological_enrich.png", width = 1390, height = 900, res = 200)
enrichplot::dotplot(enrich_GO)
dev.off()

#MF enrichment 
enrich_GO.2 <- clusterProfiler::enrichGO(gene = gene_set$ENTREZID,
                                       keyType = "ENTREZID",
                                       ont = "MF",
                                       pvalueCutoff = 0.01,
                                       pAdjustMethod = "fdr",
                                       OrgDb = "org.Hs.eg.db",
                                       universe = backgrund_set$ENTREZID
)
#save as csv file
enrich_GO.2_df <- as.data.frame(enrich_GO.2)
write.csv(enrich_GO.2_df, file = "result/Geneset_enrichment/enrich_GO_MolecularProcess_results.csv", row.names = FALSE)
print("Successfully saved Molecular process gene enrichment to result folder")

#data visualization
png("molecular_process.png", width = 1390, height = 850, res = 200)
enrichplot::dotplot(enrich_GO.2)
dev.off()

#Cellular process enrichment 
enrich_GO.3 <- enrichGO(gene = gene_set$ENTREZID,
                        ont = "CC",
                        universe = backgrund_set$ENTREZID,
                        pAdjustMethod = "fdr",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENTREZID")
#save as csv file 
enrich_GO.3_df <- as.data.frame(enrich_GO.3)
write.csv(enrich_GO.3_df, file = "result/Geneset_enrichment/enrich_GO_CellularProcess_results.csv", row.names = FALSE)
print("Successfully saved Cellular process gene enrichment to result folder")

png("cellular_process.png", width = 1390, height = 950, res = 200)
enrichplot::dotplot(enrich_GO.2)
dev.off()

cnetplot(enrich_GO, showCategory = 5)
