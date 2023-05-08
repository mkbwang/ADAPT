rm(list=ls())
library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
source('POLDA/POLDA.R')

zeller_data <- readRDS('real_data/metaanalysis/crc_zeller_results/crc_zeller_phyobj.rds')


# filter out taxa whose prevalence is smaller than 0.05
zeller_filtered_data <- filter_taxa(zeller_data, 
                                    function(x) mean(x != 0) > 0.10 , TRUE)
zeller_metadata <- sample_data(zeller_filtered_data)
zeller_count_mat <- otu_table(zeller_filtered_data)@.Data
zeller_taxonomy_table <- tax_table(zeller_filtered_data)
taxa_with_g <- zeller_taxonomy_table[, 6] !="g__"


# filter out OTU with no names on the genus level
zeller_filtered_data <- subset_taxa(zeller_filtered_data, taxa_with_g)
prevalence <- rowMeans(otu_table(zeller_filtered_data) != 0)
zeller_taxonomy_table <- tax_table(zeller_filtered_data) |> as.matrix()
zeller_metadata <- as(sample_data(zeller_filtered_data), "data.frame")
zeller_metadata$crc <- zeller_metadata$DiseaseState != "H"
sample_data(zeller_filtered_data) <- zeller_metadata
zeller_count_mat <- otu_table(zeller_filtered_data)@.Data


# Fit POLDA
start <- proc.time()
polda_result <- polda(otu_table = zeller_count_mat, metadata=zeller_metadata,
                      covar="crc", covartype="categorical",
                      startdrop="median")
polda_duration <- proc.time() - start
zeller_pval_df <- polda_result$P_Value |> na.omit()
polda_DA_taxa <- polda_result$DA_taxa
polda_reference_taxa <- polda_result$Reference_Taxa
refcounts <- zeller_count_mat[polda_reference_taxa, ]


# Fit ANCOMBC
start <- proc.time()
ancombc_result <- ancombc(zeller_filtered_data, formula = 'crc', p_adj_method='BH',
                          group='crc', struc_zero=TRUE,
                          zero_cut=0.95)
ancombc_duration <- proc.time() - start
ancombc_call <- ancombc_result$res |> as.data.frame()
colnames(ancombc_call) <- c("beta", "se", "W", "p_val", "q_val", "diff_abn")
ancombc_DA_taxa <- rownames(ancombc_call)[ancombc_call$diff_abn]



# Fit Aldex2
start <- proc.time()
aldex_result <- aldex(zeller_count_mat, zeller_metadata$crc, test='t')
aldex_duration <- proc.time() - start
aldex_DA_taxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]

zeller_taxonomy_df <- as(zeller_taxonomy_table, 'matrix') |> as.data.frame()
zeller_taxonomy_df$polda_DA <- FALSE
zeller_taxonomy_df$ancombc_DA <- FALSE
zeller_taxonomy_df$aldex_DA <- FALSE

zeller_taxonomy_df[polda_DA_taxa, 'polda_DA'] <- TRUE
zeller_taxonomy_df[ancombc_DA_taxa, 'ancombc_DA'] <- TRUE
zeller_taxonomy_df[aldex_DA_taxa, 'aldex_DA'] <- TRUE

write.csv(zeller_taxonomy_df, 'real_data/metaanalysis/zeller_performance.csv')


