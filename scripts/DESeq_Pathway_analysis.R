# ================================================================
# BRCA DIFFERENTIAL EXPRESSION & PATHWAY ANALYSIS 
# ================================================================

# ---- 1. Load Libraries ----
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)

# ---- 2. Run DESeq2 and Annotation ----
dds_corr <- DESeq(dds_corr)
results <- results(dds_corr, pAdjustMethod = "fdr", alpha = 0.05)

# Strip version numbers from Ensembl IDs for mapping (ENSG000.1 -> ENSG000)
clean_keys <- gsub("\\..*", "", rownames(results))

# Add annotation columns 
results$symbol      <- mapIds(org.Hs.eg.db, 
                              keys = clean_keys, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL", 
                              multiVals = "first")

results$ENTREZ      <- mapIds(org.Hs.eg.db, 
                              keys = clean_keys, 
                              column = "ENTREZID", 
                              keytype = "ENSEMBL", 
                              multiVals = "first")

results$description <- mapIds(org.Hs.eg.db, 
                              keys = clean_keys, 
                              column = "GENENAME", 
                              keytype = "ENSEMBL", 
                              multiVals = "first")

# Create standard significant results 
sig_results <- subset(results, padj < 0.05)

# Save progress
save(dds_corr, results, sig_results, file = "dds_results.RData")

# ---- 3. Visualization: PCA ----
dds_corr_vst <- vst(dds_corr, blind = FALSE)
plotPCA(dds_corr_vst, intgroup = "tcga.cgc_sample_sample_type", ntop = 500) +
  theme_bw() +
  geom_point(size = 3) +
  ggtitle("PCA: Top 500 variable genes")

# ---- 4. Heatmap: Healthy vs. Primary Tumor ----

# Clean group names
colData(dds_corr_vst)$tcga.cgc_sample_sample_type <- gsub(" ", "", colData(dds_corr_vst)$tcga.cgc_sample_sample_type)

# Filter and sort samples
meta_clean <- as.data.frame(colData(dds_corr_vst)) %>%
  filter(tcga.cgc_sample_sample_type != "Metastatic") %>%
  arrange(tcga.cgc_sample_sample_type) 

# Reorder matrix and take top 40 genes
mat <- assay(dds_corr_vst[row.names(sig_results), rownames(meta_clean)])[1:40, ]

# Map IDs
new_names <- sig_results[rownames(mat), "symbol"]
new_names[is.na(new_names)] <- rownames(mat)[is.na(new_names)]
rownames(mat) <- make.unique(as.character(new_names))

# Plot Heatmap
n_normal <- sum(meta_clean$tcga.cgc_sample_sample_type == "SolidTissueNormal")
pheatmap(
  mat,
  scale = "row",
  cluster_cols = FALSE,
  annotation_col = data.frame(Group = factor(meta_clean$tcga.cgc_sample_sample_type, 
                                             levels = c("SolidTissueNormal", "PrimaryTumor")),
                              row.names = rownames(meta_clean)),
  annotation_colors = list(Group = c(PrimaryTumor = "tomato", SolidTissueNormal = "green")),
  color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(255)),
  gaps_col = n_normal,
  show_colnames = FALSE,
  main = "Top 40 DEGs: Healthy vs Tumor"
)

# ---- 5. Stringent Filtering for Pathway Analysis ----
# Using your sig_3 criteria: padj < 0.01 and |LFC| > 3
sig_3_df <- as.data.frame(sig_results) %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 3) %>%
  filter(!is.na(symbol)) %>%
  distinct(symbol, .keep_all = TRUE) # Keep only unique genes

# Extract clean Entrez IDs
genes_to_test <- as.character(sig_3_df$ENTREZ[!is.na(sig_3_df$ENTREZ)])

# ---- 6. GO, KEGG, and Reactome Pathways ----

# GO Enrichment (Biological Process)
go_res <- enrichGO(gene = genes_to_test, 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   readable = TRUE)

# KEGG Enrichment
kegg_res <- enrichKEGG(gene = genes_to_test, 
                       organism = 'hsa')

# Reactome Enrichment
react_res <- enrichPathway(gene = genes_to_test, 
                           organism = "human",
                           readable = TRUE)

# ---- 7. Plotting Results ----
dotplot(go_res, 
        showCategory = 15, 
        title = "GO: Biological Process")

barplot(go_res, 
        showCategory = 10, 
        color = "p.adjust",   # Color bars by significance
        font.size = 10,
        title = "Enriched Biological Processes in BRCA - Gene Ontology")

dotplot(kegg_res,
        showCategory = 10,
        title = "Enriched Biological Processes in BRCA - KEGG")

dotplot(react_res, 
        showCategory = 15, 
        title = "Reactome Pathways")

