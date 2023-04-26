rm(list=ls())
input_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
output_folder <- '/home/wangmk/MDAWG/POLDA/simulation/ANCOMBC'

library(dplyr)
library(dacomp)

ID <- 1
input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(input_folder, input_fname))

counts_mat <- t(simulated_data$otu_tab_sim)


# find reference taxa
result.selected.references = dacomp.select_references(
  X = counts_mat, verbose = F,
  minimal_TA=10)

# correct reference taxa
result.selected.references.cleaned = dacomp.validate_references(X = counts_mat, #Counts
                                                                Y = simulated_data$metadata$X, #Traits
                                                                ref_obj = result.selected.references, #reference checked, must be object from dacomp.select_references(...)
                                                                test = DACOMP.TEST.NAME.WILCOXON, #Test used for checking, can be same test used in dacomp.test(...)
                                                                Q_validation = 0.1, #FDR level for checking
                                                                Minimal_Counts_in_ref_threshold = 1, #reference taxa will must include at least this number of reads
                                                                Reduction_Factor = 0.1, #multiplicative factor used for lowering the threshold for the number of reads required in reference taxa at each iteration
                                                                NR_perm = 1000, #number of permutations used for testing. should be at least 1/(Q_validation/ncol(X))
                                                                Verbose = T)


# DA test for the rest of taxa
result.test = dacomp.test(X = counts_mat, #counts data
                          y = simulated_data$metadata$X, #phenotype in y argument
                          # obtained from dacomp.select_references(...):
                          ind_reference_taxa = result.selected.references,
                          test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                          verbose = F)

decision_Pvals <- result.test$p.values.test
Pvals_adjusted <- p.adjust(decision_Pvals,
                                    method='BH')
Pvals_adjusted[is.na(Pvals_adjusted)] <- 1

truth <- simulated_data$taxa
truth$Taxon <- row.names(truth)

combined_results <- cbind(truth, decision_Pvals_adjusted)
library(pROC)

roc_obj <- roc(combined_results$logfold!=0, 1-combined_results$decision_Pvals_adjusted)
roc_obj$auc

