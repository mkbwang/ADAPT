rm(list=ls())
library(dplyr)
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/real_data/Yatsunenko_data'

otu_counts <- read.table(file.path(folder, 'global_gut_otu.txt'), sep='\t',
                         header=TRUE, row.names=1)
otu_counts$taxonomy <- NULL
otu_counts_mat <- as.matrix(otu_counts)

taxonomy <- read.table(file.path(folder, 'global_gut_taxonomy.txt'), sep='\t',
                       header=TRUE, row.names=1)
taxonomy_mat <- as.matrix(taxonomy)

sample_metadata <- read.table(file.path(folder, 'global_gut_metadata.txt'),
                              sep='\t', row.names=1, header=TRUE) %>%
  select(COUNTRY, AGE, SEX, FAMILY_RELATIONSHIP)
sample_metadata$AGE <- as.numeric(sample_metadata$AGE)
sample_metadata <- sample_metadata[-239,]

library(phyloseq)

phylo_data <- phyloseq(otu_table(otu_counts_mat, taxa_are_rows=TRUE),
                       tax_table(taxonomy_mat),
                       sample_data(sample_metadata))


saveRDS(phylo_data, file.path(folder, 'gut_phyloseq.rds'))


sample_metadata %>% group_by(COUNTRY) %>% summarise(count = n())


