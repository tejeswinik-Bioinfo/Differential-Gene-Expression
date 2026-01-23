library(tidyverse)
library(sva)
library(DESeq2)


countdata <- read.csv("data/raw_counts/brca_counts_data.csv", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata/brca_metadata.csv", header = TRUE, row.names = 1)

dim(countdata)
#Data cleaning

metadata_NA <- metadata %>%
  filter(!is.na(tcga.cgc_sample_sample_type))

#Filter zeroes


# Keep only genes with at least 10 counts in 10% of samples
keep_genes <- rowSums(countdata > 10) > (ncol(countdata) * 0.1)
count_filtered <- countdata[keep_genes, ]
dim(count_filtered)

save(count_filtered,countdata_clean, file = "coundata_clean_brca.RData")
# countdata_filtered_4 <- countdata[rowSums(countdata) > 124, ]
# dim(countdata_filtered_4)

countdata_clean <- count_filtered[, rownames(metadata_NA)]
identical(colnames(countdata_clean), rownames(metadata_NA))


count_matrix <- as.matrix(countdata_clean)

batch_vec <- metadata_NA$tcga.cgc_case_batch_number

condition_vec <- metadata_NA$tcga.cgc_sample_sample_type 



# We wrap this in a timing function so you know how long it took
start_time <- Sys.time()

adjusted_counts <- ComBat_seq(
  counts = count_matrix, 
  batch = metadata_NA$tcga.cgc_case_batch_number, 
  group = metadata_NA$tcga.cgc_sample_sample_type
)

end_time <- Sys.time()


dds_corr <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                                   colData = metadata_NA,
                                   design = ~ tcga.cgc_sample_sample_type)

vsd_corr <- vst(dds_corr, blind = TRUE)
plotPCA(vsd_corr, intgroup = "tcga.cgc_case_batch_number") + 
  ggtitle("PCA after ComBat-seq (Colored by Batch)")



save(adjusted_counts, file = "combat.RData")


dds_corr <- DESeq(dds_corr)






















library(DESeq2)
library(sva)
library(ggplot2)
library(dplyr)

# 1. Load Data
countdata <- read.csv("data/raw_counts/brca_counts_data.csv", header = TRUE, row.names = 1, check.names = FALSE)
metadata  <- read.csv("metadata/brca_metadata.csv", header = TRUE, row.names = 1)

# 2. Clean Metadata (Remove NAs in the condition of interest)
metadata_clean <- metadata %>%
  filter(!is.na(tcga.cgc_sample_sample_type))

# 3. Filter Genes (Keep genes with >10 counts in >10% of samples)
keep_genes <- rowSums(countdata > 10) > (ncol(countdata) * 0.1)
count_filtered <- countdata[keep_genes, ]

# 4. Final Alignment (Match columns to rows)
# Ensure only samples present in both are kept and in the same order
common_samples <- intersect(colnames(count_filtered), rownames(metadata_clean))
count_matrix   <- as.matrix(count_filtered[, common_samples])
metadata_final <- metadata_clean[common_samples, ]

stopifnot(identical(colnames(count_matrix), rownames(metadata_final)))

