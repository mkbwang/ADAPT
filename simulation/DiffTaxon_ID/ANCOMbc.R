library(ANCOMBC)
library(phyloseq)
library(dplyr)

rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'


ID <- 1

# load simulated data
simulated_data <- readRDS(file.path(folder, 'data',
                                    'semiparametric', sprintf('simulated_data_%d.rds', ID)))
# count matrix
count_mat <- otu_table(simulated_data$otu.tab.sim, taxa_are_rows = TRUE)
# binary covariate
covariate_df <- simulated_data$covariate %>% as.data.frame()
colnames(covariate_df) <- 'X'
sample_df <- sample_data(covariate_df)
test_phyloseq <- phyloseq(count_mat, sample_df)

# ancombc
ancombc_result <- ancombc(test_phyloseq, formula = 'X', p_adj_method='BH',
                          group='X', struc_zero=TRUE)


decision <- ancombc_result$res$diff_abn
colnames(decision) <- "DA"
decision$truth <- simulated_data$diff.otu.ind

write.csv(decision, file.path(folder, 'ANCOMBC_result',
                              sprintf('ANCOMBC_result_%d.csv', ID)))
