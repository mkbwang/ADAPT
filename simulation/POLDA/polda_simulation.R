rm(list=ls())

input_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data/semiparametric'

source("POLDA/POLDA.R")

ID <- 1
input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(input_folder, input_fname))

count_mat <- simulated_data$otu_tab_sim
covariate <- simulated_data$metadata
truth <- simulated_data$taxa

start <- proc.time()
pairglm_result <- pairwise_GLM(count_mat, covariate, "X")
end <- proc.time()

reftaxa <- backward_selection(count_mat, pairglm_result, start="median")

result <- reference_GLM(count_mat, covariate, "X", reftaxa)
#truth_subset <- truth %>% filter(!Taxa %in% reftaxa)
truth$Taxon <- rownames(truth)
combined_result <- truth %>% left_join(result, by="Taxon")
combined_result$pval_adjust[is.na(combined_result$pval_adjust)] <- 1
table(combined_result$logfold!=0, combined_result$pval_adjust<0.05)

library(pROC)

roc_obj <- roc(combined_result$logfold != 0, 1-combined_result$pval_adjust)
roc_obj$auc



