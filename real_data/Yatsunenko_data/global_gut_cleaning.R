# clean data and retain counts on the phylum level

rm(list=ls())

library(readxl)
library(tidyverse)
library(nloptr)
library(phyloseq)
library(stringr)
library(ggpubr)
library(magrittr)


# metadata read in

meta_data=read_tsv("real_data/Yatsunenko_data/global_gut_metadata.txt")
meta_data=meta_data%>%transmute(Sample.ID=`SampleID`, age=AGE, sex=SEX, country=COUNTRY)%>%
  arrange(Sample.ID) %>% as.data.frame()
meta_data=meta_data[complete.cases(meta_data), ]
meta_data$age=as.numeric(meta_data$age)
meta_data$country=recode(meta_data$country, `GAZ:Malawi` = "MA",
                         `GAZ:United States of America` = "US", `GAZ:Venezuela` = "VEN")
rownames(meta_data) <- meta_data$Sample.ID
meta_data$Sample.ID <- NULL

# metadata read in 2

meta_data2=read_excel("real_data/Yatsunenko_data/NIHMS365354-supplement-3.xls", sheet = 1, skip = 2)
meta_data2=meta_data2[!is.na(meta_data2$`Sample Identifier`), ]
meta_data2=meta_data2%>%
  transmute(Sample.ID=`Sample Identifier`, age=`Age (Years)`, gender=Gender,
            bmi=`BMI (kg/m2)`, breast.fed=`Breast-fed`, country=Country)%>%
  arrange(Sample.ID)
meta_data2$age=signif(as.numeric(meta_data2$age), digits = 2)
meta_data2$bmi=signif(as.numeric(meta_data2$bmi), digits = 2)
meta_data2$breast.fed[which(meta_data2$breast.fed=="NA1")]="NA"

# read in taxonomy
taxonomy=read_tsv("real_data/Yatsunenko_data/global_gut_taxonomy.txt") |> as.data.frame()
rownames(taxonomy) <- as.character(taxonomy$OTU_ID)
taxonomy$OTU_ID <- NULL


# read in the OTU table and aggregate into phylum level
otu_table=read_tsv("real_data/Yatsunenko_data/global_gut_otu.txt") |> as.data.frame()
rownames(otu_table) <- as.character(otu_table$OTU_ID)
otu_table=otu_table[, -c(1, 532)] # get rid of taxonomy column


# select a subset of taxa and individuals as template
US_adult <- meta_data[meta_data$country == "US" & meta_data$age > 18, ] |> na.omit()
otu_table_subset <- otu_table[, rownames(US_adult)]

taxa_meancount <- rowMeans(otu_table_subset)
taxa_prevalence <- rowMeans(otu_table_subset != 0)
taxa_filter <- taxa_prevalence > 0.5 & taxa_meancount > 10

otu_table_subset <- otu_table_subset[taxa_filter, ]
subset_taxonomy <- taxonomy[taxa_filter, ]
simulation_template <- phyloseq(otu_table(as.matrix(otu_table_subset), taxa_are_rows=TRUE),
                                tax_table(as.matrix(subset_taxonomy)),
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

