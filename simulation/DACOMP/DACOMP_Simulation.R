rm(list=ls())
input_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data/semiparametric'
output_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/DACOMP_result/Simulation'

library(dplyr)
library(dacomp)

ID <- 1
input_fname <- sprintf('simulated_data_%d.RDS', ID)
simulated_data <- readRDS(file.path(input_folder, input_fname))

counts_mat <- t(simulated_data$otu_tab_sim)

# find reference taxa
result.selected.references = dacomp.select_references(
  X = counts_mat, verbose = F)

# DA test for the rest of taxa
result.test = dacomp.test(X = counts_mat, #counts data
                          y = simulated_data$metadata$X, #phenotype in y argument
                          # obtained from dacomp.select_references(...):
                          ind_reference_taxa = result.selected.references,
                          test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                          verbose = F)

decision_Pvals <- result.test$p.values.test
decision_Pvals[is.na(decision_Pvals)] <- 1

truth <- simulated_data$taxa
subset_truth <- truth$logfold[!is.na(decision_Pvals)]

library(pROC)

roc_obj <- roc(subset_truth != 0, 1-decision_Pvals[!is.na(decision_Pvals)])
roc_obj$auc

