rm(list=ls())
library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
source('POLDA/POLDA.R')
baxter_data <- readRDS('real_data/metaanalysis/crc_baxter_results/crc_baxter_phyobj.rds')


# filter out taxa whose prevalence is smaller than 0.05
baxter_filtered_data <- filter_taxa(baxter_data, 
                                      function(x) mean(x != 0) > 0.05 , TRUE)
baxter_metadata <- sample_data(baxter_filtered_data)
prevalence <- rowMeans(otu_table(baxter_filtered_data) != 0)
baxter_count_mat <- otu_table(baxter_filtered_data)@.Data
baxter_taxonomy_table <- tax_table(baxter_filtered_data)
taxa_with_g <- baxter_taxonomy_table[, 6] !="g__"

# filter out OTU with no names on the genus level
baxter_filtered_data <- subset_taxa(baxter_filtered_data, taxa_with_g)
prevalence <- rowMeans(otu_table(baxter_filtered_data) != 0)
baxter_taxonomy_table <- tax_table(baxter_filtered_data) |> as.matrix()
baxter_metadata <- as(sample_data(baxter_filtered_data), "data.frame")
baxter_metadata$crc <- baxter_metadata$DiseaseState != "H"
sample_data(baxter_filtered_data) <- baxter_metadata
baxter_count_mat <- otu_table(baxter_filtered_data)@.Data



# Fit POLDA
start <- proc.time()
polda_result <- polda(otu_table = baxter_count_mat, metadata=baxter_metadata,
                      covar="crc", covartype="categorical",
                      startdrop="median")
polda_DA_taxa <- polda_result$DA_taxa
polda_duration <- proc.time() - start



# Fit ANCOMBC
start <- proc.time()
ancombc_result <- ancombc(baxter_filtered_data, formula = 'crc', p_adj_method='BH',
                          group='crc', struc_zero=TRUE,
                          zero_cut=0.95)
ancombc_duration <- proc.time() - start
ancombc_call <- ancombc_result$res$diff_abn
ancombc_DA_taxa <- rownames(ancombc_call)[ancombc_call$crcTRUE]


# Fit Aldex2
start <- proc.time()
aldex_result <- aldex(baxter_count_mat, baxter_metadata$crc, test='t')
aldex_duration <- proc.time() - start

aldex_DA_taxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]

baxter_taxonomy_df <- as(baxter_taxonomy_table, 'matrix') |> as.data.frame()
baxter_taxonomy_df$polda_DA <- FALSE
baxter_taxonomy_df$ancombc_DA <- FALSE
baxter_taxonomy_df$aldex_DA <- FALSE

baxter_taxonomy_df[polda_DA_taxa, 'polda_DA'] <- TRUE
baxter_taxonomy_df[ancombc_DA_taxa, 'ancombc_DA'] <- TRUE
baxter_taxonomy_df[aldex_DA_taxa, 'aldex_DA'] <- TRUE

write.csv(baxter_taxonomy_df, 'real_data/metaanalysis/baxter_performance.csv')

