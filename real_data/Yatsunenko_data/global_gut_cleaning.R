# clean data and retain counts on the phylum level

rm(list=ls())

library(readxl)
library(tidyverse)
library(nloptr)
library(phyloseq)
library(stringr)
library(ggpubr)
library(magrittr)


# metadata read in 2

meta_data=read_excel("real_data/Yatsunenko_data/NIHMS365354-supplement-3.xls", sheet = 1, skip = 2)
meta_data=meta_data[!is.na(meta_data$`Sample Identifier`), ]
meta_data=meta_data%>%
  transmute(Sample.ID=`Sample Identifier`, age=`Age (Years)`, gender=Gender,country=Country,
            Depth = `Number of V4-16S rRNA Illumina sequences`)%>%
  arrange(Sample.ID)
meta_data$age=signif(as.numeric(meta_data$age), digits = 2)


# read in taxonomy
# taxonomy=read_tsv("real_data/Yatsunenko_data/global_gut_taxonomy.txt") |> as.data.frame()
# rownames(taxonomy) <- as.character(taxonomy$OTU_ID)
# taxonomy$OTU_ID <- NULL


# read in the OTU table and aggregate into phylum level
otu_table=read_tsv("real_data/Yatsunenko_data/Yatsunenko2012.txt") |> as.data.frame()
rownames(otu_table) <- as.character(otu_table$OTU_ID)
taxonomy <- otu_table[, 530]
otu_table=otu_table[, -c(1, 530)] # get rid of taxonomy column


# select a subset of taxa and individuals as template
## the subset of individuals and US adults
US_adult <- meta_data[meta_data$country == "USA" & meta_data$age > 18, ] |> na.omit() |>
  as.data.frame()
rownames(US_adult) <- US_adult$Sample.ID
US_adult$Sample.ID <- NULL
otu_table_subset <- otu_table[, rownames(US_adult)] |> as.matrix()
taxonomies <- sapply(taxonomy, function(fullname) strsplit(fullname, split=";")[[1]]) |>
  t()
rownames(taxonomies) <- rownames(otu_table)
colnames(taxonomies) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

## the subset of taxa are those whose median relative abundance ranks top 1000
depths <- colSums(otu_table_subset)
depths_mat <- t(replicate(nrow(otu_table_subset), depths))
relative_abundance_mat <- otu_table_subset / depths_mat
median_relative_abundance <- apply(relative_abundance_mat, 1, median)
relabd_rank <- length(median_relative_abundance) + 1 - rank(median_relative_abundance)
selected_taxa <- which(relabd_rank < 1000)


otu_table_subset <- otu_table_subset[selected_taxa, ]
subset_taxonomies <- taxonomies[selected_taxa, ]
simulation_template <- phyloseq(otu_table(as.matrix(otu_table_subset), taxa_are_rows=TRUE),
                                tax_table(as.matrix(subset_taxonomies)),
                                sample_data(US_adult))
saveRDS(simulation_template, "simulation/data/Yatsunenko_template.rds")

otu_table$OTU_ID=taxonomy$Phylum[match(otu_table$OTU_ID, taxonomy$OTU_ID)]
phylum_table=otu_table%>%group_by(OTU_ID)%>%
  summarise_all(sum)
non_info_pos=grep("\\p__\\b", phylum_table$OTU_ID) # Exact match
phylum_table=phylum_table[-non_info_pos, ]
phylum_table=as.data.frame(phylum_table)


rownames(phylum_table) <- phylum_table$OTU_ID
phylum_table <- phylum_table[, -1]

meta_data <- as.data.frame(meta_data)
# there is an error in metadata
meta_data$country[47] <- "VEN"
row.names(meta_data) <- meta_data$Sample.ID
meta_data$Sample.ID <- NULL

write.csv(phylum_table,
          "real_data/Yatsunenko_data/phylum_table.csv")
write.csv(meta_data, "real_data/Yatsunenko_data/cleaned_metadata.csv")


data_subset <- function(otu_count, metadata, select_country, agegroup="infant"){

  if (agegroup == "infant"){
    metadata_subset <- metadata%>%filter(country%in%select_country, age<=2)
  } else{
    metadata_subset <- metadata%>%filter(country%in%select_country, age>=18, age<=40)
  }

  otu_count_subset <- otu_count[, row.names(metadata_subset)]
  prevalence <- rowMeans(otu_count_subset != 0)
  selected_taxa <- names(prevalence[prevalence > 0.2])
  otu_count_subset <- otu_count_subset[selected_taxa, ]

  phy_obj <- phyloseq(otu_table(otu_count_subset, taxa_are_rows=TRUE),
                      sample_data(metadata_subset))
  list(metadata_subset = metadata_subset,
       otu_count_subset = otu_count_subset,
       phy_obj = phy_obj)

}

# US vs Venezuela
US_VEN_infant <- data_subset(phylum_table, meta_data, c("US", "VEN"), "infant")
write.csv(US_VEN_infant$metadata_subset, 'real_data/Yatsunenko_data/US_Venezuela/US_VEN_infant_metadata.csv')
write.csv(US_VEN_infant$otu_count_subset, 'real_data/Yatsunenko_data/US_Venezuela/US_VEN_infant_OTU.csv')
saveRDS(US_VEN_infant$phy_obj, file='real_data/Yatsunenko_data/US_Venezuela/US_VEN_infant_phyloseq.rds')


US_VEN_adult <- data_subset(phylum_table, meta_data, c("US", "VEN"), "adult")
write.csv(US_VEN_adult$metadata_subset, 'real_data/Yatsunenko_data/US_Venezuela/US_VEN_adult_metadata.csv')
write.csv(US_VEN_adult$otu_count_subset, 'real_data/Yatsunenko_data/US_Venezuela/US_VEN_adult_OTU.csv')
saveRDS(US_VEN_adult$phy_obj, file='real_data/Yatsunenko_data/US_Venezuela/US_VEN_adult_phyloseq.rds')


# US vs Malawi
US_MA_infant <- data_subset(phylum_table, meta_data, c("US", "MA"), "infant")
write.csv(US_MA_infant$metadata_subset, 'real_data/Yatsunenko_data/US_Malawi/US_MA_infant_metadata.csv')
write.csv(US_MA_infant$otu_count_subset, 'real_data/Yatsunenko_data/US_Malawi/US_MA_infant_OTU.csv')
saveRDS(US_MA_infant$phy_obj, file='real_data/Yatsunenko_data/US_Malawi/US_MA_infant_phyloseq.rds')

US_MA_adult <- data_subset(phylum_table, meta_data, c("US", "MA"), "adult")
write.csv(US_MA_adult$metadata_subset, 'real_data/Yatsunenko_data/US_Malawi/US_MA_adult_metadata.csv')
write.csv(US_MA_adult$otu_count_subset, 'real_data/Yatsunenko_data/US_Malawi/US_MA_adult_OTU.csv')
saveRDS(US_MA_adult$phy_obj, file='real_data/Yatsunenko_data/US_Malawi/US_MA_adult_phyloseq.rds')


# Malawi vs Venezuela
MA_VEN_infant <- data_subset(phylum_table, meta_data, c("VEN", "MA"), "infant")
write.csv(MA_VEN_infant$metadata_subset, 'real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_infant_metadata.csv')
write.csv(MA_VEN_infant$otu_count_subset, 'real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_infant_OTU.csv')
saveRDS(MA_VEN_infant$phy_obj, file='real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_infant_phyloseq.rds')

MA_VEN_adult <- data_subset(phylum_table, meta_data, c("VEN", "MA"), "adult")
write.csv(MA_VEN_adult$metadata_subset, 'real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_adult_metadata.csv')
write.csv(MA_VEN_adult$otu_count_subset, 'real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_adult_OTU.csv')
saveRDS(MA_VEN_adult$phy_obj, file='real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_adult_phyloseq.rds')

