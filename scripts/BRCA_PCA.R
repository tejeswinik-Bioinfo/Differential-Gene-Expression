

BiocManager::install("recount3")
library(recount3)

########################################################
###Download dataset using recount3
########################################################

human_projects <- recount3::available_projects()
proj_info <- subset(human_projects, project == "BRCA" & project_type == "data_sources")
rse_brca <- recount3::create_rse(proj_info)


head(rse_brca)
class(rse_brca)
#[1] "RangedSummarizedExperiment"
#attr(,"package")
#[1] "SummarizedExperiment"
typeof(rse_brca)
#[1] "S4"

#since it is a s4 object , cannot write to .csv. S4 object contains the sample info, gene info, metadata.
#Extract the metadata and sample,genes separatey and write it to .csv


########################################################
###Extracting count data & metadata and writing to .csv
########################################################

#Extracting the count matrix from rse_brca
brca_count_data <- assay(rse_brca)
write.csv(as.data.frame(brca_count_data), "data/brca_counts_data.csv", row.names = TRUE)
#Extracting the metadata from rse_brca
brca_metadata <- colData(rse_brca)
write.csv(as.data.frame(brca_metadata),"brca_metadata.csv", row.names = TRUE) 

# #Assigning a column name
# colnames(brca_metadata)[1] <- "sample_ID"

#Extracting required column from the metadata - creating a subset of a metadata


keep_cols <- c("external_id",
               "tcga.cgc_sample_sample_type",
               "tcga.gdc_platform",
               "tcga.gdc_cases.demographic.gender",
               "tcga.gdc_cases.demographic.race",
               "tcga.gdc_cases.demographic.ethnicity",
               "tcga.gdc_cases.tissue_source_site.project",
               "tcga.gdc_cases.tissue_source_site.name",
               "tcga.gdc_cases.diagnoses.tumor_stage",
               "tcga.gdc_cases.samples.sample_id",
               "tcga.gdc_cases.samples.is_ffpe",
               "tcga.gdc_cases.samples.days_to_collection",
               "tcga.cgc_sample_country_of_sample_procurement",
               "tcga.cgc_case_year_of_diagnosis",
               "tcga.cgc_case_age_at_diagnosis",
               "tcga.cgc_case_batch_number")

clean_metadata <- as.data.frame(colData(rse_brca)[, keep_cols])
colnames(clean_metadata) <- c("ID", "Sample_type", "Platform", "Gender", "Race", "Ethnicity", "tissue_source_site", "source_site_name", "tumor_stage",
                              "sample_ID", "is_ffpe", "days_to_collection", "country_of_sample_procurement", "year_of_diagnosis", "age_at_diagnosis","batch_number")
write.csv(clean_metadata, "metadata/subset_brca_metadata.csv", row.names = TRUE)

###################################################################################################
#Load the data
counts <- read.csv("data/raw_counts/brca_counts_data.csv", header =TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata/subset_brca_metadata.csv", header = TRUE, row.names = 1)

# if(ncol(counts) == nrow(metadata)) {
#   print("Sample counts match")
# }
# 
# identical(rownames(metadata), colnames(counts))

#Matching the rownames of metadata and the column names of count data



# 
# all(rownames(metadata) %in% colnames(counts))
# counts <- counts[, rownames(metadata)]

library(tidyverse)
metadata_NA <- metadata %>%
  filter(!is.na(Sample_type))

counts_NA <- counts[ , rownames(metadata_NA)]
metadata <- metadata[match(colnames(counts_NA), metadata_NA$ID), ]
identical(rownames(metadata_NA), colnames(counts_NA))

########################################################
###Creating DESeq object, normalization & PCA plot
########################################################


library(DESeq2)

#creating dds object

dds <- DESeqDataSetFromMatrix(countData = counts_NA,
                              colData = metadata_NA,
                              design = ~Sample_type)


#Normalization

brca_vsd <- vst(dds, blind = TRUE)


library(ggplot2)

#PCA plot - Before batch effect correction and after normalization

plotPCA(brca_vsd, intgroup = "Sample_type") +
  geom_point(size = 4) +
  theme_light() +
  coord_cartesian(ylim = c(-50,50)) +
  labs(
    title = "PRINCIPAL COMPONENT ANALYSIS",
    subtitle = "Samples clustered by Gene Expression",
    caption = "Before Batch effect correction - Data Normalized"
    )

#Write Column names to .csv
metadata_colnames <- colnames(brca_metadata) #column names of original metadata
subset_metadata_colnames <- colnames(metadata_NA) #column names of subset metadata

write.csv(metadata_colnames, "metadata/brca_metadata_colnames.csv", row.names = TRUE)
write.csv(subset_metadata_colnames, "metadata/brca_metadata_subset_colnames.csv", row.names = TRUE)

#PCA plot - After batch effect correction and after normalization
