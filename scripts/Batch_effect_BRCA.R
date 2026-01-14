metadata <- read.csv("metadata/brca_metadata.csv", header = TRUE, row.names = 1)

# Assuming your metadata is called 'metadata'
# 1. Create a summary table of all columns
df_summary <- data.frame(
  column_name = colnames(metadata),
  unique_values = sapply(metadata, function(x) length(unique(x))),
  na_count = sapply(metadata, function(x) sum(is.na(x))),
  data_type = sapply(metadata, class)
)
write.csv(df_summary, "Results/summary_brca.csv", row.names = TRUE)

# 2. Filter for "Batch Candidates"
# A good batch candidate usually has > 1 unique value but < 20 (to avoid IDs)
# and is not mostly NAs.
batch_candidates <- df_summary[df_summary$unique_values > 1 & 
                                 df_summary$unique_values < 10 & 
                                 df_summary$na_count < (nrow(metadata) * 0.5), ]

print(batch_candidates)

# 3. Search for specific keywords (like 'date', 'center', 'site', 'platform')
grep("date|center|site|platform|batch|year", colnames(metadata), value = TRUE, ignore.case = TRUE)

###############################################

subset_metadata <- read.csv("metadata/subset_brca_metadata.csv", header = TRUE, row.names = 1)
# 1. See all unique locations and their sample counts
table(subset_metadata$source_site_name)

# 2. To see it as a percentage (helpful for finding dominant sites)
prop.table(table(subset_metadata$source_site_name)) * 100

# 3. To see ALL features for ALL columns in your subset
# This will print the unique counts for every column in your list
lapply(subset_metadata[, c("source_site_name", "Platform", "is_ffpe", "batch_number")], table)





library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)

# Step 1: Prepare data (Assuming 'counts' is your matrix and 'metadata' is your df)
# We use VST to stabilize variance for PCA
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~)
vsd <- vst(dds, blind = TRUE)
pca_data <- prcomp(t(assay(vsd)))

# Step 2: Extract top 5 Principal Components
pcs <- pca_data$x[, 1:5]

# Step 3: Correlate PCs with Metadata
# We use Kruskal-Wallis for categorical variables to get a p-value of association
target_metadata <- metadata[, c(bio_var, batch_vars)]

calc_association <- function(pc_vector, meta_vector) {
  if(is.numeric(meta_vector)) {
    return(cor(pc_vector, meta_vector, method = "spearman"))
  } else {
    # Return 1 - p-value so that high values = strong association
    test <- kruskal.test(pc_vector ~ as.factor(meta_vector))
    return(1 - test$p.value) 
  }
}

# Build the association matrix
assoc_matrix <- apply(pcs, 2, function(pc) {
  apply(target_metadata, 2, function(meta) calc_association(pc, meta))
})

# Step 4: Plot the Heatmap
Heatmap(assoc_matrix, 
        name = "1 - p-value", 
        column_title = "Principal Components",
        row_title = "Metadata Variables",
        col = colorRampPalette(c("white", "red"))(100),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", assoc_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
        })




countdata <- read.csv("data/raw_counts/brca_counts_data.csv" , header =TRUE, row.names = 1)

rownames(metadata) <- metadata$external_id

# Reorder metadata to match count matrix
metadata <- metadata[match(colnames(countdata), metadata$external_id), ]


#PCA plot
library(DESeq2)
library(ggplot2)
library(patchwork) # For combining plots

# 1. Prepare the DESeq2 object and transform
# Use the raw counts you got from recount3
dds <- DESeqDataSetFromMatrix(countData = countdata, 
                              colData = metadata, 
                              design = ~ tcga.cgc_sample_sample_type)

colnames(countdata)
metadata$external_id
rownames(metadata)



# VST is faster than rlog for large datasets like TCGA (1200+ samples)
vsd <- vst(dds, blind = TRUE)

# 2. Extract PCA data
# We use intgroup to tell DESeq2 which columns to keep for the plot
pcaData_bio <- plotPCA(vsd, intgroup = "tcga.cgc_sample_sample_type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData_bio, "percentVar"))

# 3. Create Plot 1: Colored by Biology (Sample Type)
p1 <- ggplot(pcaData_bio, aes(PC1, PC2, color = tcga.cgc_sample_sample_type)) +
  geom_point(size = 3, alpha = 0.7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Colored by Condition (Biology)") +
  theme_minimal()

# 4. Create Plot 2: Colored by Suspected Batch (Batch Number)
# Note: Ensure the batch column name matches your metadata exactly
pcaData_batch <- plotPCA(vsd, intgroup = "tcga.cgc_case_batch_number", returnData = TRUE)

p2 <- ggplot(pcaData_batch, aes(PC1, PC2, color = as.factor(tcga.cgc_case_batch_number))) +
  geom_point(size = 3, alpha = 0.7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Colored by TCGA Batch Number") +
  theme_minimal() +
  theme(legend.position = "none") # Hide legend if there are too many batches (38+)

# 5. Display side-by-side
p1 + p2




dim(countdata)
dim(metadata)














