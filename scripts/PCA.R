#Load Libraries

library(ggplot2)
library(DESeq2)

#Load count data file
data <- read.csv("GSE159337_rare_dis_read_count_ori.csv", header = TRUE, row.names = 1)

#Load metadata
metadata <- read.delim("human_rna_seq/hyperlipidemia_18.12.2025/metadata.txt", header = TRUE)

#convert Sample_ID to rownames
rownames(metadata)<- metadata$Sample_Name

#Match the read count data and metadata
metadata<- metadata[match(colnames(data), metadata$Sample_Name), ]
head(metadata)

#Creating DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~isolate)

#Normalization
dds_rlog <- rlog(dds, blind =TRUE)

# #PCA plot before normalization
# pca_raw <- ggplot(dds, aes(x = PC1, y = PC2, color = isolate)) +
#   theme_minimal() +
#   theme_light() +
#   geom_point(size = 4) +
#   labs(
#     title = "Principal Component Analysis (PCA)",
#     subtitle = "Gene Expression Profiles",
#     caption = "Raw counts"
#   )


#PCA plot
plotPCA(dds_rlog, intgroup = "isolate")+
  geom_point(size = 4) +
  theme_light() +
  #theme_minimal() +
  coord_cartesian(ylim = c(-20, 20)) +
  scale_color_manual(values = c("healthy_control" = "blue", "patient" = "orange")) +
  labs(
      title = "Principal Component Analysis (PCA)",
      subtitle = "Samples clustered by gene Expression Profiles",
      caption = "Data normalized using rlog"
      )