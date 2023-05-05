rm(list=ls())

library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
source('POLDA/POLDA.R')
edd_singh_data <- readRDS('real_data/metaanalysis/edd_singh_results/edd_singh_phyobj.rds')

# filter out taxa whose prevalence is smaller than 5%
edd_singh_filtered_data <- filter_taxa(edd_singh_data, function(x) mean(x != 0) > 0.05 , TRUE)
edd_singh_count_mat <- otu_table(edd_singh_filtered_data)@.Data
edd_singh_taxonomy_table <- tax_table(edd_singh_filtered_data)
taxa_with_g <- edd_singh_taxonomy_table[, 6] !="g__"


# filter out OTU with no names on the genus level
edd_singh_filtered_data <- subset_taxa(edd_singh_filtered_data, taxa_with_g)
edd_singh_taxonomy_table <- tax_table(edd_singh_filtered_data)
edd_singh_metadata <- as(sample_data(edd_singh_filtered_data), "data.frame")
edd_singh_metadata$diarrhea <- edd_singh_metadata$DiseaseState != "H"
sample_data(edd_singh_filtered_data) <- edd_singh_metadata
edd_singh_count_mat <- otu_table(edd_singh_filtered_data)@.Data


# Fit POLDA
start <- proc.time()
polda_result <- polda(otu_table = edd_singh_count_mat, metadata=edd_singh_metadata,
                      covar="diarrhea", covartype="categorical",
                      startdrop="median")
polda_duration <- proc.time() - start
polda_DA_taxa <- c(polda_result$DA_taxa, polda_result$Structural_zero_Taxa)


refcounts <- edd_singh_count_mat[polda_result$Reference_Taxa, ]


# Fit ANCOMBC
start <- proc.time()
ancombc_result <- ancombc(edd_singh_filtered_data, formula = 'diarrhea', p_adj_method='BH',
                          group='diarrhea', struc_zero=TRUE,
                          zero_cut=0.95)
ancombc_duration <- proc.time() - start
ancombc_call <- ancombc_result$res$diff_abn
ancombc_DA_taxa <- rownames(ancombc_call)[ancombc_call$diarrheaTRUE]


# Fit Aldex2
start <- proc.time()
aldex_result <- aldex(edd_singh_count_mat, edd_singh_metadata$diarrhea, test='t')
aldex_duration <- proc.time() - start
aldex_DA_taxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]


edd_singh_taxonomy_df <- as(edd_singh_taxonomy_table, 'matrix') |> as.data.frame()
edd_singh_taxonomy_df$polda_DA <- FALSE
edd_singh_taxonomy_df$ancombc_DA <- FALSE
edd_singh_taxonomy_df$aldex_DA <- FALSE

edd_singh_taxonomy_df[polda_DA_taxa, 'polda_DA'] <- TRUE
edd_singh_taxonomy_df[ancombc_DA_taxa, 'ancombc_DA'] <- TRUE
edd_singh_taxonomy_df[aldex_DA_taxa, 'aldex_DA'] <- TRUE

write.csv(edd_singh_taxonomy_df, 'real_data/metaanalysis/edd_singh_performance.csv')

