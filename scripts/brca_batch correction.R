
###########################################################################
##Step 1: Data Preparation & Filtering
###########################################################################

# 1. Load Data
countdata <- read.csv("data/raw_counts/brca_counts_data.csv", header = TRUE, row.names = 1, check.names = FALSE)
metadata  <- read.csv("metadata/brca_metadata.csv", header = TRUE, row.names = 1)

# 2. Clean Metadata (Remove NAs in the condition of interest)
metadata_clean <- metadata %>%
  filter(!is.na(tcga.cgc_sample_sample_type))

# 3. Filter Genes (Keep genes with >10 counts in >10% of samples)
keep_genes <- rowSums(countdata > 10) > (ncol(countdata) * 0.1)
count_filtered <- countdata[keep_genes, ]

# 4. Match columns to rows

common_samples <- intersect(colnames(count_filtered), rownames(metadata_clean))
count_matrix   <- as.matrix(count_filtered[, common_samples])
metadata_final <- metadata_clean[common_samples, ]

stopifnot(identical(colnames(count_matrix), rownames(metadata_final)))

###########################################################################
##Step 2: Before Batch Correction
###########################################################################

# Create DESeq object for Raw Data
dds_raw <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData   = metadata_final,
                                  design    = ~ tcga.cgc_sample_sample_type)

# Variance Stabilizing Transformation (VST)
vsd_raw <- vst(dds_raw, blind = TRUE)

# Prepare Plot Data
pca_raw_df <- plotPCA(vsd_raw, intgroup = "tcga.cgc_case_batch_number", returnData = TRUE)
percentVar <- round(100 * attr(pca_raw_df, "percentVar"))

# Plot 1: Before Correction
p1 <- ggplot(pca_raw_df, aes(PC1, PC2, color = as.factor(tcga.cgc_case_batch_number))) +
  geom_point(size = 2, alpha = 0.6) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("1. BEFORE ComBat-seq (Colored by Batch)") +
  theme_minimal() + theme(legend.position = "none")

###########################################################################
##Step 3: After Batch Correction
###########################################################################

start_t <- Sys.time()

# Run the progress-aware version of the function you edited
adjusted_counts <- ComBat_seq(
  counts = count_matrix, 
  batch  = metadata_final$tcga.cgc_case_batch_number, 
  group  = metadata_final$tcga.cgc_sample_sample_type
)

end_t <- Sys.time()


# Save adjusted counts immediately to avoid re-running
save(adjusted_counts, file = "brca_adjusted_counts.RData")

# --- Visualizing Adjusted Data ---
dds_corr <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                                   colData   = metadata_final,
                                   design    = ~ tcga.cgc_sample_sample_type)

vsd_corr <- vst(dds_corr, blind = TRUE)
pca_corr_df <- plotPCA(vsd_corr, intgroup = "tcga.cgc_case_batch_number", returnData = TRUE)
percentVar_corr <- round(100 * attr(pca_corr_df, "percentVar"))

# Plot 2: After Correction
p2 <- ggplot(pca_corr_df, aes(PC1, PC2, color = as.factor(tcga.cgc_case_batch_number))) +
  geom_point(size = 2, alpha = 0.6) +
  xlab(paste0("PC1: ", percentVar_corr[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_corr[2], "% variance")) +
  ggtitle("2. AFTER ComBat-seq (Colored by Batch)") +
  theme_minimal() + theme(legend.position = "none")

# Display plots side-by-side
p1 + p2

