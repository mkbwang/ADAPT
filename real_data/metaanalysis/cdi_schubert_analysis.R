rm(list=ls())
library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
source('POLDA/POLDA.R')
schubert_data <- readRDS('real_data/metaanalysis/cdi_schubert_results/cdi_schubert_phyobj.rds')

# filter out taxa whose prevalence is smaller than 0.05
schubert_filtered_data <- filter_taxa(schubert_data, 
                                      function(x) mean(x != 0) > 0.05 , TRUE)
schubert_count_mat <- otu_table(schubert_filtered_data)@.Data
schubert_taxonomy_table <- tax_table(schubert_filtered_data)
taxa_with_g <- schubert_taxonomy_table[, 6] !="g__"

# filter out OTU with no names on the genus level
schubert_filtered_data <- subset_taxa(schubert_filtered_data, taxa_with_g)
schubert_taxonomy_table <- tax_table(schubert_filtered_data)
schubert_metadata <- as(sample_data(schubert_filtered_data), "data.frame")
schubert_metadata$diarrhea <- schubert_metadata$DiseaseState != "H"
sample_data(schubert_filtered_data) <- schubert_metadata
schubert_count_mat <- otu_table(schubert_filtered_data)@.Data

# Fit POLDA
start <- proc.time()
polda_result <- polda(otu_table = schubert_count_mat, metadata=schubert_metadata,
                      covar="diarrhea", covartype="categorical",
                      startdrop="median")
polda_duration <- proc.time() - start
polda_DA_taxa <- c(polda_result$DA_taxa, polda_result$Structural_zero_Taxa)

tau_density_fit <- density(polda_result$Tau_hat,
                           bw="SJ", kernel="epanechnikov")
tau_mode <- tau_density_fit$x[which.max(tau_density_fit$y)]

# Fit ANCOMBC
start <- proc.time()
ancombc_result <- ancombc(schubert_filtered_data, formula = 'diarrhea', p_adj_method='BH',
                          group='diarrhea', struc_zero=TRUE,
                          zero_cut=0.95)
ancombc_duration <- proc.time() - start
ancombc_call <- ancombc_result$res$diff_abn
ancombc_DA_taxa <- rownames(ancombc_call)[ancombc_call$diarrheaTRUE]

# Fit Aldex2
start <- proc.time()
aldex_result <- aldex(schubert_count_mat, schubert_metadata$diarrhea, test='t')
aldex_duration <- proc.time() - start

aldex_DA_taxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]

# locom gets stuck, need to double check
# locom_result <- locom(otu.table=t(schubert_count_mat), 
#                       Y=as.numeric(schubert_metadata$diarrhea),
#                       fdr.nominal=0.2,
#                       prev.cut=0.05,
#                       adjustment="BH",
#                       n.cores = 4)

schubert_taxonomy_df <- as(schubert_taxonomy_table, 'matrix') |> as.data.frame()
schubert_taxonomy_df$polda_DA <- FALSE
schubert_taxonomy_df$ancombc_DA <- FALSE
schubert_taxonomy_df$aldex_DA <- FALSE

schubert_taxonomy_df[polda_DA_taxa, 'polda_DA'] <- TRUE
schubert_taxonomy_df[ancombc_DA_taxa, 'ancombc_DA'] <- TRUE
schubert_taxonomy_df[aldex_DA_taxa, 'aldex_DA'] <- TRUE

write.csv(schubert_taxonomy_df, 'real_data/metaanalysis/schubert_performance.csv')
