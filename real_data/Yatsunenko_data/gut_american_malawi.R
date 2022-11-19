rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/real_data/Yatsunenko_data'
# folder = '.'
library(dplyr)
library(phyloseq)

data <- readRDS(file.path(folder, 'gut_phyloseq.rds'))
otu_mat <- otu_table(data)
taxa_info <- tax_table(data)
indv_metadata <- sample_data(data)

# select american and malawi infants

indv_IDs <- rownames(indv_metadata)
age_filter <- indv_metadata$AGE < 2
country_filter <- indv_metadata$COUNTRY %in% c("GAZ:Malawi", "GAZ:United States of America")

selected_individuals <- indv_IDs[age_filter & country_filter]
metadata <- data.frame(ID = selected_individuals,
                       Country = indv_metadata[selected_individuals, "COUNTRY"]) %>%
  na.omit()
selected_individuals <- metadata$ID
metadata$indicator <- metadata$COUNTRY == "GAZ:Malawi"

americans_id <- selected_individuals[!metadata$indicator]
malawi_id <- selected_individuals[metadata$indicator]


otu_subset <- otu_mat[, selected_individuals] %>% as.matrix()

# total_prevalence <- apply(otu_subset,1,function(x){
#   subp <- mean(x > 0)
# })
#
#
# taxa_mask <- total_prevalence >= 0.25
#
# prevalences <- data.frame(american = american_prevalence[taxa_mask],
#                           malawi = malawi_prevalence[taxa_mask])
#
# counts <- otu_subset[taxa_mask, ] %>% as.matrix()


source('/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/ancom_simple.R')

preprocess_result <- feature_table_pre_process(otu_subset, metadata,
                                               'ID', 'indicator', out_cut = 0.05,
                                               zero_cut=0.75, lib_cut=1e5, neg_lb = TRUE)

filtered_counts <- preprocess_result$feature_table
# don't agree with the outlier zeros
filtered_counts[is.na(filtered_counts)] <- 0

filtered_taxa <- row.names(filtered_counts)

metadata <- preprocess_result$meta_data
struc_zeros <- preprocess_result$structure_zeros

# taxa that have been identified to be differentialy abundant
diff_taxa <- filtered_taxa[which((struc_zeros[, 1] == 1) | (struc_zeros[, 2] == 1))]
retained_taxa <- filtered_taxa[which((struc_zeros[, 1] == 0) & (struc_zeros[, 2] == 0))]

final_filtered_counts <- filtered_counts[retained_taxa, ]

cleaned_gut_data <- list(struc_zero_taxa = diff_taxa,
                         counts = final_filtered_counts,
                         metadata = metadata)

saveRDS(cleaned_gut_data, file.path(folder, 'subset_gut.rds'))

