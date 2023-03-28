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
  arrange(Sample.ID)
meta_data=meta_data[complete.cases(meta_data), ]
meta_data$age=as.numeric(meta_data$age)
meta_data$country=recode(meta_data$country, `GAZ:Malawi` = "MA",
                         `GAZ:United States of America` = "US", `GAZ:Venezuela` = "VEN")

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

write.csv(meta_data2, "real_data/Yatsunenko_data/cleaned_metadata.csv", row.names=FALSE)

# read in taxonomy
taxonomy=read_tsv("real_data/Yatsunenko_data/global_gut_taxonomy.txt")
taxonomy=taxonomy%>%rowwise()%>%
  mutate(genus_name=paste(Phylum, Genus, sep = ";"))


# read in the OTU table and aggregate into phylum level
otu_table=read_tsv("real_data/Yatsunenko_data/global_gut_otu.txt")
otu_table=otu_table[, -532] # get rid of taxonomy column
otu_table$OTU_ID=taxonomy$Phylum[match(otu_table$OTU_ID, taxonomy$OTU_ID)]
phylum_table=otu_table%>%group_by(OTU_ID)%>%
  summarise_all(sum)
non_info_pos=grep("\\p__\\b", phylum_table$OTU_ID) # Exact match
phylum_table=phylum_table[-non_info_pos, ]
phylum_table=as.data.frame(phylum_table)

write.csv(phylum_table,
          "real_data/Yatsunenko_data/phylum_table.csv",
          row.names=FALSE)

# summarize metadata for two age groups: 0-2 and 18-40
meta_data_infant_US_VEN = meta_data2%>%filter(country%in%c("US", "VEN"), age<=2)
# meta_data_18_40 = meta_data2%>%filter(age>=18, age<=40)


obs.abn <- phylum_table
rownames(obs.abn)=obs.abn$OTU_ID
obs.abn=obs.abn[, -1]

prevalence <- rowMeans(obs.abn == 0)
selected_taxa <- names(prevalence[prevalence < 0.9])

# load functions for ANCOMBC
source('ANCOM_ANCOMBC/ancombc.R')



