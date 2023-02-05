
library(tidyverse)
library(blockmodels)
library(ggplot2)
library(knitr)
theme_set(theme_bw())

rm(list=ls())


source('simulation/DiffTaxon_ID/SBM_utils.R')


# for (replicate in 1:100){
  # load the truth
replicate <- 10
data_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data/semiparametric'
simulated_data = readRDS(file.path(data_folder,
                                   sprintf("simulated_data_%d.rds", replicate)))

taxa_info <- data.frame(Taxa = simulated_data$otu.names,
                              DA = simulated_data$diff.otu.ind)

# load the GLM pairwise analysis result
GLM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/glmdisp_result/semiparametric'
GLM_result <- read.csv(file.path(GLM_folder, 'taxapair',
                                 sprintf('glmdisp_result_%d.csv', replicate)))
# GLM_result$adjusted_pval <- p.adjust(GLM_result$pval, method='BH')
GLM_result$decision <- GLM_result$pval < 0.05
GLM_decision <- GLM_result %>% dplyr::select(T1, T2, decision)
GLM_decision_wide <- GLM_decision %>% reshape(idvar="T1", timevar="T2", direction="wide")
dimension <- nrow(GLM_decision_wide) + 1

rm(GLM_result)
rm(GLM_decision)

## change from long vector to matrix
GLM_decision_mat <- matrix(FALSE, nrow=dimension,
                           ncol=dimension)

for (j in 1:(dimension-1)){
  GLM_decision_mat[j, (j+1):dimension] <- GLM_decision_wide[j, (j+1):dimension] %>%
    as.numeric()
}
GLM_decision_mat <- GLM_decision_mat + t(GLM_decision_mat)
# plotMyMatrix(GLM_decision_mat, dimLabels =c('Taxa'))


# fit SBM model
init_membership_degree <- cap_membership(degree_init(GLM_decision_mat))
init_membership_sc <- spectral_clustering(GLM_decision_mat)
init_membership_flat <- flat_init(GLM_decision_mat)
result_degree <- complete_EM(GLM_decision_mat, init_membership_degree)
result_sc <- complete_EM(GLM_decision_mat, init_membership_sc)
result_flat <- complete_EM(GLM_decision_mat, init_membership_flat)
result_pred_prob_degree <- result_degree$membership
result_pred_prob_sc <- result_sc$membership

# simplest SBM model with binary entries
taxaSBM <- BM_bernoulli(membership_type="SBM_sym", adj=GLM_decision_mat,
                        explore_min=2, explore_max=2, plotting='', verbosity=6)
taxaSBM$estimate()

result_pred_prob_package <- taxaSBM$memberships[[2]]$Z
difference <- result_pred_prob_sc - result_pred_prob_package
max(abs(difference))

## make the decision
taxa_info$prediction <- result_degree$membership[, 1] > 0.5

output_filename <- sprintf('glmdisp_SBM_decision_%d.csv', replicate)
write.csv(taxa_foldchange, file.path(GLM_folder, 'taxonindv', output_filename),
          row.names = FALSE)
# }


# for (replicate in 1:100){
#   print(replicate)
#   # load the truth
#   data_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data'
#   simulated_data = readRDS(file.path(data_folder,
#                                      sprintf("simulated_data_%d.rds", replicate)))
#   abn_info <- simulated_data$mean.eco.abn
#   taxa_foldchange <- data.frame(Taxa = row.names(abn_info),
#                                 effect.size = abn_info$effect.size)
#
#   # load the GLM pairwise analysis result
#   ancom_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/ancom_result'
#   ancom_result = read.csv(file.path(ancom_folder,
#                                   sprintf("ancom_result_%d.csv", replicate)))
#   ancom_result$adjusted_pval <- p.adjust(ancom_result$pval, method='BH')
#   ancom_result$decision <- ancom_result$adjusted_pval < 0.05
#   ancom_decision <- ancom_result %>% dplyr::select(T1, T2, decision)
#   ancom_decision_wide <- ancom_decision %>% reshape(idvar="T1", timevar="T2", direction="wide")
#   dimension <- nrow(ancom_decision_wide) + 1
#
#   rm(ancom_result)
#   rm(ancom_decision)
#
#   ## change from long vector to matrix
#   ancom_decision_mat <- matrix(FALSE, nrow=dimension,
#                              ncol=dimension)
#
#   for (j in 1:(dimension-1)){
#     ancom_decision_mat[j, (j+1):dimension] <- ancom_decision_wide[j, (j+1):dimension] %>%
#       as.numeric()
#   }
#   ancom_decision_mat <- ancom_decision_mat + t(ancom_decision_mat)
#   # plotMyMatrix(GLM_decision_mat, dimLabels =c('Taxa'))
#
#   # simplest SBM model with binary entries
#   taxaSBM <- ancom_decision_mat %>%
#     estimateSimpleSBM("bernoulli", dimLabels = 'Taxa',
#                       estimOptions = list(verbosity = 1, plot = FALSE, nbCores=4,
#                                           exploreMin=2, exploreMax=4, nbBlocksRange=c(2, 4)))
#
#
#   taxaSBM$setModel(2)
#
#   ## make the decision by checking the connection expectation within own group
#   ## low expectation corresponds to nonDA group
#   expectations <- diag(taxaSBM$expectation)
#   taxa_foldchange$prediction <- expectations == max(expectations)
#   taxa_foldchange$truth <- taxa_foldchange$effect.size !=1
#
#   output_filename <- sprintf('ancom_SBM_decision_%d.csv', replicate)
#   write.csv(taxa_foldchange, file.path(ancom_folder, 'SBM', output_filename),
#             row.names = FALSE)
# }
#
#
