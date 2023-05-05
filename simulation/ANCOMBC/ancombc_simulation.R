library(ANCOMBC)
library(phyloseq)
library(dplyr)

rm(list=ls())
input_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
output_folder <- '/home/wangmk/MDAWG/POLDA/simulation/ANCOMBC'

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# load simulated data
simulated_data <- readRDS(file.path(input_folder, sprintf('simulated_data_%d.rds', ID)))
# count matrix
count_mat <- otu_table(simulated_data$otu_tab_sim, taxa_are_rows = TRUE)
# binary covariate
covariate_df <- simulated_data$metadata
sample_df <- sample_data(covariate_df)
test_phyloseq <- phyloseq(count_mat, sample_df)

# ancombc
ancombc_result <- ancombc(test_phyloseq, formula = 'X', p_adj_method='BH',
                          group='X', struc_zero=FALSE)


Q_vals <- ancombc_result$res$q_val
Q_vals$Taxon <- row.names(Q_vals)
truth <- simulated_data$taxa
truth$Taxon <- row.names(truth)

combined_result <- truth %>% left_join(Q_vals, by='Taxon')
colnames(combined_result)[3] <- "Pval_Adjust"
combined_result$Taxon <- rownames(combined_result)

write.csv(combined_result, 
          file.path(output_folder, sprintf('ancombc_simulation_result_%d.csv', ID)),
          row.names=FALSE)

