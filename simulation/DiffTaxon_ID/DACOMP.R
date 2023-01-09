rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'

library(dplyr)
library(dacomp)

ID <- 1
input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(folder, 'data', 'semiparametric', input_fname))

counts_mat <- t(simulated_data$otu.tab.sim)

# find reference taxa
result.selected.references = dacomp.select_references(
  X = counts_mat, verbose = F)

# DA test for the rest of taxa
result.test = dacomp.test(X = counts_mat, #counts data
                          y = as.vector(simulated_data$covariate), #phenotype in y argument
                          # obtained from dacomp.select_references(...):
                          ind_reference_taxa = result.selected.references,
                          test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                          verbose = F)

decision <- result.test$dsfdr_rejected

result_df <- data.frame(DA=FALSE, truth=simulated_data$diff.otu.ind)
result_df$DA[decision] <- TRUE
row.names(result_df) <- colnames(counts_mat)

output_fname <- sprintf('ALDEx2_%d.csv', ID)
write.csv(result_df, file.path(folder, 'ALDEx2_result', output_fname))

