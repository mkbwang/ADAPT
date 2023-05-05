rm(list=ls())

input_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
output_folder <- '/home/wangmk/MDAWG/POLDA/simulation/ANCOM'

source("/home/wangmk/MDAWG/POLDA/simulation/ANCOM/ancombc_prep.R")
source("/home/wangmk/MDAWG/POLDA/simulation/ANCOM/ancom_utils.R")

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# load simulated data
simulated_data <- readRDS(file.path(input_folder, sprintf('simulated_data_%d.rds', ID)))
# count matrix
count_mat <- otu_table(simulated_data$otu_tab_sim, taxa_are_rows = TRUE)
# binary covariate
covariate_df <- simulated_data$metadata
covariate_df$X <- as.factor(covariate_df$X)
sample_df <- sample_data(covariate_df)
test_phyloseq <- phyloseq(count_mat, sample_df)

out <- ancom(data = test_phyloseq, assay_name = "counts",
             tax_level = NULL, phyloseq = NULL,
             p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
             main_var = "X", adj_formula = NULL,
             rand_formula = NULL, lme_control = NULL,
             struc_zero = TRUE, neg_lb = FALSE, alpha = 0.05, n_cl = 1)

result <- out$res 
truth <- simulated_data$taxa
truth$taxon <- row.names(truth)

performance <- truth %>% left_join(result, by="taxon")

write.csv(performance, 
        file.path(output_folder, sprintf('ancom_simulation_result_%d.csv', ID)),
        row.names=FALSE)

