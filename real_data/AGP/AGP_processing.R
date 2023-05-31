# Simulation Template Generation
# The purpose of this script is to generate a OTU table template for simulation studies
# The template is based on a recently uploaded batch of American Gut Project 16S data (Oct 2022)
# Qiita project ID 10317, sample ID 13742

rm(list=ls())
library(biomformat)
library(phyloseq)

# load downloded data from Qiita
AGP_batch <- read_biom('real_data/AGP/AGP_13742_otu.biom')
count_mat <- biom_data(AGP_batch) |> as.matrix()
taxonomy_info <- observation_metadata(AGP_batch)
metadata_info <- sample_metadata(AGP_batch)

# Ckeck out the sequencing depth distribution
depths <- colSums(count_mat)
hist(depths, nclass=20)

# I notice a bimodal pattern, so I only keep the samples whose sequencing depths are larger than 50000
selected_samples <- names(depths[depths > 5e4])
num_sample <- length(selected_samples)
subset_count_mat <- count_mat[, selected_samples]

# retain taxa with top 1000 average prevalence(existence across all samples)
taxa_prevalence <- rowMeans(subset_count_mat != 0)
sorted_prevalence <- sort(taxa_prevalence, decreasing=TRUE)
selected_taxa <- names(sorted_prevalence[1:1000])
num_taxa <- length(selected_taxa)

subset_count_mat <- count_mat[selected_taxa, selected_samples]
subset_depths <- colSums(subset_count_mat)
hist(subset_depths, nclass=20)
subset_taxonomy <- taxonomy_info[selected_taxa,]


# observe the sequencing depths and relative abundances
relative_abundance <- subset_count_mat / t(replicate(num_taxa, subset_depths))
max_rel_abd <- apply(relative_abundance, 1, max)
med_rel_abd <- apply(relative_abundance, 1, median)

# wrap into phyloseq and save
colnames(subset_count_mat) <- substring(colnames(subset_count_mat), 7)
subset_metadata <- data.frame(X=c(rep(0, (num_sample-1)/2), rep(1, (num_sample+1)/2)))
rownames(subset_metadata) <- colnames(subset_count_mat)

template_phyobj <- phyloseq(otu_table(subset_count_mat, taxa_are_rows=TRUE),
                            tax_table(as.matrix(subset_taxonomy)),
                            sample_data(subset_metadata))

saveRDS(template_phyobj, 'simulation/data/AGP_template.rds')

