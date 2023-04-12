library(dada2)
rm(list=ls())
memory.limit(25000) # scale memory
library(biomformat)
library(phyloseq)
library(gtools)

file_path <- "real_data/Crohn_Vandeputte/gut_otu_table.RData" #load counts data
load(file_path)
# Average_Cell_Count = read.csv(file = 'real_data/Crohn_Vandeputte/FlowCytometryCounts.csv')$Average_Cell_Count #load flow cytometry counts
Y = rep(0,95) #group labeling is such that the first 29 subjects are CD, the rest are healthy.
Y[1:29] = 1

#show a box plot for the distribution of average cell counts across subjects
#boxplot(Average_Cell_Count[Y==1],Average_Cell_Count[Y==0],names = c('Sick','Healthy'),main='Boxplot for cell counts in Healthy/Sick')
#mean(Average_Cell_Count[Y==0])/(mean(Average_Cell_Count[Y==1])) #ratio between mean number of reads - discussed in paper
#2.508608
#median(Average_Cell_Count[Y==0])/(median(Average_Cell_Count[Y==1])) # ratio between median number of reads
#3.086826

# order OTU table (from deblur) by original sample names
reordering_permutation = order(colnames(otu_table))
otu_table = otu_table[, reordering_permutation]
colnames(otu_table) #it is now ordered


set.seed(1)
seqs <- rownames(otu_table)
taxa <- assignTaxonomy(seqs,
                       "real_data/Crohn_Vandeputte/gg_13_8_train_set_97.fa.gz",
                       multithread=T,
                       verbose = T,
                       minBoot = 80)


otu_df <-  as.data.frame(otu_table)

otu_df$genus <- as.vector(taxa[, "Genus"])
otu_df$genus <- gsub("\\[|\\]", "", otu_df$genus)
library(dplyr)
otu_df_subset <- otu_df %>% na.omit() %>% filter(genus != "g__")

genus_counts <- otu_df_subset %>% group_by(genus) %>% summarise(count=n())
otu_df_genus_level <- otu_df_subset %>% group_by(genus) %>%
  summarise_all(sum) %>% as.data.frame()

#prevalence of taxa across subjects
rownames(otu_df_genus_level) <- otu_df_genus_level$genus
otu_df_genus_level$genus <- NULL
otu_genux_mat <- as.matrix(otu_df_genus_level)
prevalence_matrix = 1*(otu_genux_mat>0)
prev_taxa <- as.numeric(apply(prevalence_matrix,1,mean))
hist(prev_taxa)


crohn_group <- otu_genux_mat[, seq(1,29)]
crohn_existence <- 1*(crohn_group > 0)
crohn_prev <- as.numeric(apply(crohn_existence,1,sum))

control_group <- otu_genux_mat[, seq(30, 95)]
control_existence <- 1*(control_group > 0)
control_prev <- as.numeric(apply(control_existence,1,sum))

valid <- (crohn_prev > 2) & (control_prev > 4)

selected_genus <- rownames(otu_df_genus_level)[valid]

otu_genus_subset <- otu_df_genus_level[valid, ]
# taxonomy_subset <- taxa[valid, ]
metadata <- data.frame(Crohn = Y)
rownames(metadata) <- colnames(otu_genus_subset)


# permute the counts of each taxon to generate null cases
perm_func <- function(otutable, group1=29, group2=66){

  otutable_permuted <- otutable
  for (j in 1:nrow(otutable)){

    # find all nonzero values
    allvals <- otutable[j, ]
    nonzero_vals <- allvals[allvals != 0]
    num_nonzero <- length(nonzero_vals)

    # allocation of nonzero values in the two groups
    g1alloc <- rbinom(1, size=num_nonzero, prob=group1/(group1+group2))
    if (g1alloc < max(2, num_nonzero - group2)){
      g1alloc <- max(2, num_nonzero - group2)
    } else if (g1alloc > min(group1, num_nonzero - 2)){
      g1alloc <- min(group1, num_nonzero - 2)
    }

    g2alloc <- num_nonzero - g1alloc

    # maintain the number of nonzero entries in each group
    g1_indices <- sample(seq(1, group1), g1alloc)
    g2_indices <- sample(seq(group1+1, group1+group2), g2alloc)

    nonzero_indices <- c(g1_indices, g2_indices)
    otutable_permuted[j, ] <- 0
    otutable_permuted[j, nonzero_indices] <- permute(nonzero_vals)
  }
  return(otutable_permuted)
}


otu_genus_subset_perm1 <- perm_func(otu_genus_subset)
otu_genus_subset_perm2 <- perm_func(otu_genus_subset)
otu_genus_subset_perm3 <- perm_func(otu_genus_subset)



# crohn_group <- otu_table_subset_perm3[, seq(1,29)]
# crohn_existence <- 1*(crohn_group > 0)
# crohn_prev <- as.numeric(apply(crohn_existence,1,sum))
#
#
# control_group <- otu_table_subset_perm3[, seq(30, 95)]
# control_existence <- 1*(control_group > 0)
# control_prev <- as.numeric(apply(control_existence,1,sum))
#
# # rename the taxa
# rownames(otu_table_subset) <- sprintf("OTU%d", seq(1, nrow(otu_table_subset)))
# rownames(otu_table_subset_perm1) <- sprintf("OTU%d", seq(1, nrow(otu_table_subset)))
# rownames(otu_table_subset_perm2) <- sprintf("OTU%d", seq(1, nrow(otu_table_subset)))
# rownames(otu_table_subset_perm3) <- sprintf("OTU%d", seq(1, nrow(otu_table_subset)))
# rownames(taxonomy_subset) <- sprintf("OTU%d", seq(1, nrow(otu_table_subset)))

phy_obj <- phyloseq(otu_table(otu_genus_subset, taxa_are_rows=TRUE),
                    sample_data(metadata))
phy_permuted1 <- phyloseq(otu_table(otu_genus_subset_perm1, taxa_are_rows=TRUE),
                          sample_data(metadata))
phy_permuted2 <- phyloseq(otu_table(otu_genus_subset_perm2, taxa_are_rows=TRUE),
                          sample_data(metadata))
phy_permuted3 <- phyloseq(otu_table(otu_genus_subset_perm3, taxa_are_rows=TRUE),
                          sample_data(metadata))
saveRDS(phy_obj, "real_data/Crohn_Vandeputte/Crohn_phyloseq.rds")
saveRDS(phy_permuted1, "real_data/Crohn_Vandeputte/Crohn_phyloseq_permuted1.rds")
saveRDS(phy_permuted2, "real_data/Crohn_Vandeputte/Crohn_phyloseq_permuted2.rds")
saveRDS(phy_permuted3, "real_data/Crohn_Vandeputte/Crohn_phyloseq_permuted3.rds")


# phy_obj <- phyloseq(otu_table(otu_df_genus_level_subset, taxa_are_rows=TRUE),
#                     sample_data(metadata))


# otu_df <- as.data.frame(otu_table_subset)
# otu_perm1_df <- as.data.frame(otu_table_subset_perm1)
# otu_perm2_df <- as.data.frame(otu_table_subset_perm2)
# otu_perm3_df <- as.data.frame(otu_table_subset_perm3)

write.csv(otu_genus_subset, "real_data/Crohn_Vandeputte/Crohn_otu.csv")
write.csv(otu_genus_subset_perm1, "real_data/Crohn_Vandeputte/Crohn_otu_perm1.csv")
write.csv(otu_genus_subset_perm2, "real_data/Crohn_Vandeputte/Crohn_otu_perm2.csv")
write.csv(otu_genus_subset_perm3, "real_data/Crohn_Vandeputte/Crohn_otu_perm3.csv")
write.csv(metadata, "real_data/Crohn_Vandeputte/Crohn_metadata.csv")
