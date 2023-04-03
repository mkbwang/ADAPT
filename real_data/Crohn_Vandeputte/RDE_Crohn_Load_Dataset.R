library(dada2)
rm(list=ls())
memory.limit(25000) # scale memory
library(biomformat)
library(phyloseq)
Q_LEVEL = 0.1 #BH threshold level

file_path <- "real_data/Crohn_Vandeputte/gut_otu_table.RData" #load counts data
load(file_path)
Average_Cell_Count = read.csv(file = 'real_data/Crohn_Vandeputte/FlowCytometryCounts.csv')$Average_Cell_Count #load flow cytometry counts
Y = rep(0,95) #group labeling is such that the first 29 subjects are CD, the rest are healthy.
Y[1:29] = 1

#show a box plot for the distribution of average cell counts across subjects
boxplot(Average_Cell_Count[Y==1],Average_Cell_Count[Y==0],names = c('Sick','Healthy'),main='Boxplot for cell counts in Healthy/Sick')
mean(Average_Cell_Count[Y==0])/(mean(Average_Cell_Count[Y==1])) #ratio between mean number of reads - discussed in paper
#2.508608
median(Average_Cell_Count[Y==0])/(median(Average_Cell_Count[Y==1])) # ratio between median number of reads
#3.086826

# order OTU table (from deblur) by original sample names
reordering_permutation = order(colnames(otu_table))
otu_table = otu_table[, reordering_permutation]
colnames(otu_table) #it is now ordered

set.seed(1)
otu_df <-  as.data.frame(otu_table)
seqs <- rownames(otu_df)
taxa <- assignTaxonomy(seqs,
                       "real_data/Crohn_Vandeputte/gg_13_8_train_set_97.fa.gz",
                       multithread=T,
                       verbose = T,
                       minBoot = 80)

otu_df$genus <- as.vector(taxa[, "Genus"])
otu_df$genus <- gsub("\\[|\\]", "", otu_df$genus)
library(dplyr)
otu_df_subset <- otu_df %>% na.omit() %>% filter(genus != "g__")

genus_counts <- otu_df_subset %>% group_by(genus) %>% summarise(count=n())

otu_df_genus_level <- otu_df_subset %>% group_by(genus) %>% summarise_all(sum) %>% as.data.frame()
row.names(otu_df_genus_level) <- otu_df_genus_level$genus
otu_df_genus_level$genus <- NULL

prevalence <- rowMeans(otu_df_genus_level != 0)

otu_df_genus_level_subset <- otu_df_genus_level[prevalence > 0.2, ]

write.csv(otu_df_genus_level_subset, "real_data/Crohn_Vandeputte/Crohn_otu.csv")

metadata <- data.frame(Crohn = Y)
rownames(metadata) <- colnames(otu_df_genus_level_subset)

write.csv(metadata, "real_data/Crohn_Vandeputte/Crohn_metadata.csv")

phy_obj <- phyloseq(otu_table(otu_df_genus_level_subset, taxa_are_rows=TRUE),
                    sample_data(metadata))

saveRDS(phy_obj, "real_data/Crohn_Vandeputte/Crohn_phyloseq.rds")


