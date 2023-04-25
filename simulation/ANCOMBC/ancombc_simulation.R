library(ANCOMBC)
library(phyloseq)
library(dplyr)

rm(list=ls())
input_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data/semiparametric'
output_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/ANCOMBC_result/Simulation'


ID <- 1

# load simulated data
simulated_data <- readRDS(file.path(input_folder, sprintf('simulated_data_%d.RDS', ID)))
# count matrix
count_mat <- otu_table(simulated_data$otu_tab_sim, taxa_are_rows = TRUE)
# binary covariate
covariate_df <- simulated_data$metadata
sample_df <- sample_data(covariate_df)
test_phyloseq <- phyloseq(count_mat, sample_df)

# ancombc
ancombc_result <- ancombc(test_phyloseq, formula = 'X', p_adj_method='BH',
                          group='X', struc_zero=TRUE)


Q_vals <- ancombc_result$res$q_val
truth <- simulated_data$taxa
combined_result <- cbind(truth, Q_vals)
table(combined_result$logfold!=0, combined_result$X<0.16)

library(pROC)

roc_obj <- roc(truth$logfold != 0, 1-Q_vals$X)
roc_obj$auc

