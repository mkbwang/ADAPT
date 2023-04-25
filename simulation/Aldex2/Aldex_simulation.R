rm(list=ls())
library(dplyr)
library(ALDEx2)

input_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data/semiparametric'
output_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/ALDEx2_result/Simulation'

ID <- 1
input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(input_folder, input_fname))

count_mat <- simulated_data$otu_tab_sim
covariate <- as.vector(simulated_data$metadata$X)
truth <- simulated_data$taxa

aldex_result <- aldex(count_mat, covariate, test='t')

result_df <- data.frame(DA_Pval=aldex_result$wi.eBH,
                        DA_truth=truth$logfold !=0)
table(result_df$DA_truth, result_df$DA_Pval < 0.05)

library(pROC)

roc_obj <- roc(result_df$DA_truth, 1-result_df$DA_Pval)
roc_obj$auc

