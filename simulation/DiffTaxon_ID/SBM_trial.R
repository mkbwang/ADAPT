
library(tidyverse)
library(blockmodels)
library(ggplot2)
library(knitr)
theme_set(theme_bw())

rm(list=ls())




n <- 12 # nodes
Z <- matrix(0, nrow=n, ncol=2)
Z[1:6, 1] <- 1; Z[7:12, 2] <- 1
# Z<-diag(Q)%x%matrix(1,npc,1)
P<-matrix(c(0.6, 0.9, 0.9, 0.2), nrow=2)
P[lower.tri(P)]<-t(P)[lower.tri(P)]
set.seed(20)
M<-1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z)) ## adjacency matrix
M[lower.tri(M)]<-t(M)[lower.tri(M)]
# M <- matrix(1, nrow=10, ncol=10)
diag(M) <- 0

## estimation
my_model <- BM_bernoulli("SBM_sym", M , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=6)
my_model$estimate()
my_model$memberships[[2]]$Z
which.max(my_model$ICL)



# example_mat <- rbind(c(0,1,0,1,0,1),
#                      c(0,0,0,1,0,0),
#                      c(0,1,0,1,0,1),
#                      c(0,0,0,0,0,1),
#                      c(0,1,0,1,0,1),
#                      c(0,0,0,0,0,0))
# example_mat <- example_mat + t(example_mat)
# # original_plot <- plotMyMatrix(example_mat, dimLabels =c('example'))
# exampleSBM <- BM_bernoulli(membership_type="SBM", adj=example_mat,
#                            explore_min=2, explore_max=2)
# exampleSBM$estimate()
#
# membershipprob <- exampleSBM$memberships[[2]]$Z

# exampleSBM$plot(type="data")
# clustered_plot <- plot(exampleSBM, type = "data", dimLabels  = c('example'))

# library(cowplot)
# toyexample_plot <- plot_grid(original_plot, clustered_plot, ncol=2)
#

# for (replicate in 1:100){
  # load the truth
replicate <- 0
data_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data/semiparametric'
simulated_data = readRDS(file.path(data_folder,
                                   'mini_simulated_data.rds'))

taxa_info <- data.frame(Taxa = simulated_data$otu.names,
                              DA = simulated_data$diff.otu.ind)

# load the GLM pairwise analysis result
GLM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/glmdisp_result/semiparametric'
GLM_result <- read.csv(file.path(GLM_folder, 'taxapair', 'glmdisp_result_mini.csv'))
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

# simplest SBM model with binary entries
taxaSBM <- BM_bernoulli(membership_type="SBM_sym", adj=GLM_decision_mat,
                        explore_min=2, explore_max=2, plotting='', verbosity=6)
taxaSBM$estimate()


## make the decision by checking the connection expectation within own group
## low expectation corresponds to nonDA group
expectations <- taxaSBM$memberships[[2]]$Z
taxa_info$prediction <- expectations[, 1] < expectations[, 2]

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
