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
zero_counts <- rowSums(baxter_count_mat == 0)

# three steps for polda
baxter_pairglm <- pairwise_GLM(baxter_count_mat,
                               baxter_metadata,
                               covar="crc")

reftaxa_search <- backward_selection(baxter_count_mat,
                                     baxter_pairglm,
                                     start="median")


polda_result <- reference_GLM(baxter_count_mat, baxter_metadata, 'crc',
                              reftaxa=reftaxa_search$reftaxa)


# Fit ANCOMBC
ancombc_result <- ancombc(baxter_filtered_data, formula = 'crc', p_adj_method='BH',
                          group='crc', struc_zero=TRUE,
                          zero_cut=0.95)
ancombc_call <- ancombc_result$res |> as.data.frame()
colnames(ancombc_call) <- c("beta", "se", "W", "p_val", "q_val", "diff_abn")
ancombc_DA_taxa <- rownames(ancombc_call)[ancombc_call$crcTRUE]


# Fit Aldex2
aldex_result <- aldex(baxter_count_mat, baxter_metadata$crc, test='t')
aldex_DA_taxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]

baxter_taxonomy_df <- as(baxter_taxonomy_table, 'matrix') |> as.data.frame()
baxter_taxonomy_df$polda_effect <- NA
baxter_taxonomy_df[polda_result$Taxon, 'polda_effect'] <- polda_result$effect
baxter_taxonomy_df$polda_rawp <- NA
baxter_taxonomy_df[polda_result$Taxon, 'polda_rawp'] <- polda_result$pval
baxter_taxonomy_df$polda_adjustedp <- NA
baxter_taxonomy_df[polda_result$Taxon, 'polda_adjustedp'] <- polda_result$pval_adjust
baxter_taxonomy_df$ancombc_rawp <- ancombc_call$p_val
baxter_taxonomy_df$ancombc_adjustedp <- ancombc_call$q_val
baxter_taxonomy_df$aldex_rawp <- aldex_result$wi.ep
baxter_taxonomy_df$aldex_adjustedp <- aldex_result$wi.eBH


write.csv(baxter_taxonomy_df, 'real_data/metaanalysis/baxter_performance.csv')

