rm(list=ls())
library(ANCOMBC)
library(phyloseq)

data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
ancom_folder <- '/home/wangmk/MDAWG/POLDA/simulation/ANCOM'

source("/home/wangmk/MDAWG/POLDA/simulation/ANCOM/ancom_utils.R")

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 10

# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata_null <- AGP_null$sample_metadata
phyobj_null <- phyloseq(otu_table(AGP_null$count_mat, taxa_are_rows=TRUE),
                        sample_data(sample_metadata_null))

out <- ancom(data = phyobj_null, assay_name = "counts",
             p_adj_method = "BH", prv_cut = 0.05, lib_cut = 0,
             main_var = "X", struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05, n_cl = 2)


# load simulated data
simulated_data <- readRDS(file.path(input_folder, sprintf('simulated_data_%d.rds', ID)))
# count matrix
count_mat <- otu_table(simulated_data$otu_tab_sim, taxa_are_rows = TRUE)
# binary covariate
covariate_df <- simulated_data$metadata
covariate_df$X <- as.factor(covariate_df$X)
sample_df <- sample_data(covariate_df)
test_phyloseq <- phyloseq(count_mat, sample_df)


result <- out$res 
truth <- simulated_data$taxa
truth$taxon <- row.names(truth)

performance <- truth %>% left_join(result, by="taxon")

write.csv(performance, 
        file.path(output_folder, sprintf('ancom_simulation_result_%d.csv', ID)),
        row.names=FALSE)

