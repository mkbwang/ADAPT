
library(tidyverse)
library(blockmodels)
library(ggplot2)
library(cowplot)
library(knitr)
theme_set(theme_bw())

rm(list=ls())


source('simulation/DiffTaxon_ID/SBM_utils.R')


# for (replicate in 1:100){
  # load the truth

replicate <- 20
data_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data/semiparametric'
GLM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/glmdisp_result/semiparametric'

GLM_result <- read.csv(file.path(GLM_folder, 'taxapair',
                                 sprintf('glmdisp_result_null_%d.csv', replicate)))


simulated_data <- readRDS(file.path(data_folder,
                                   sprintf("simulated_data_null_%d.rds", replicate)))

## set threshold at 0.05
GLM_result$link05 <- GLM_result$pval < 0.05
GLM_decision_mat_high <- adjmat_generation(data_folder, GLM_folder, id=replicate, threshold=0.05)
## set threshold at 0.01
GLM_result$link01 <- GLM_result$pval < 0.01
GLM_decision_mat_low <- adjmat_generation(data_folder, GLM_folder, id=replicate, threshold=0.01)


# visualize the histograms of degrees
degrees_high <- colMeans(GLM_decision_mat_high)
degrees_low <- colMeans(GLM_decision_mat_low)
degrees_df <- data.frame(Degrees=c(degrees_high, degrees_low),
                         Cutoff = rep(c("Cutoff=0.05", "Cutoff=0.01"), each=300))
degrees_df_nonzero <- degrees_df %>% filter(Degrees > 0)
ggplot(degrees_df_nonzero %>% filter(Cutoff == "Cutoff=0.05"), aes(x=Degrees)) +
  geom_histogram(binwidth=0.02, color="black", fill="white")+
  xlab("Degrees/300 (Only Non-zero degrees included)")+ylab("Count")


# visualize the adjacency matrix
name_order_1 <- sprintf("Taxon%d", seq(1, 300))
name_order_2 <- sprintf("Taxon%d", seq(300, 1))


GLM_result_copy_1 <- GLM_result %>% select(T1, T2, link05, link01) %>%
  rename(from=T1, to=T2) %>% mutate(to=factor(to, levels=name_order_1),
                                    from=factor(from, levels=name_order_2))


GLM_result_copy_2 <- GLM_result_copy_1
GLM_result_copy_2$from <- GLM_result_copy_1$to
GLM_result_copy_2$to <- GLM_result_copy_1$from
GLM_result_copy_3 <- GLM_result_copy_2
GLM_result_copy_3$to <- GLM_result_copy_3$from
GLM_result_copy_3$link05 <- FALSE
GLM_result_copy_3$link01 <- FALSE
GLM_result_copy <- rbind(GLM_result_copy_1, GLM_result_copy_2, GLM_result_copy_3)


ggplot(GLM_result_copy, aes(x = from, y = to, fill = link05)) +
  geom_raster() +
  theme_bw() + scale_fill_manual(values = c("white", "black"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none")

ggplot(GLM_result_copy, aes(x = from, y = to, fill = link01)) +
  geom_raster() +
  theme_bw() + scale_fill_manual(values = c("white", "black"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none")




# fit SBM model
init_membership_degree <- cap_membership(degree_init(GLM_decision_mat_high))
init_membership_sc <- spectral_clustering(GLM_decision_mat_high)
init_membership_flat <- flat_init(GLM_decision_mat_high)


result_degree <- complete_EM(GLM_decision_mat_high, init_membership_degree)
result_sc <- complete_EM(GLM_decision_mat_high, init_membership_sc)
result_flat <- complete_EM(GLM_decision_mat_high, init_membership_flat)


# plot trajectories
membership_trajectories_summary <- function(VI_result){
  num_iter <- length(VI_result$membership)
  num_taxa <- nrow(VI_result$membership[[1]])
  membership_trajectories <- data.frame(Iter = rep(1:num_iter, each=num_taxa),
                                               Taxon = rep(1:num_taxa, num_iter),
                                               Membership = 0)
  membership_trajectories$Membership <- do.call(c, lapply(VI_result$membership,
                                                                 function(element) element[, 1]))
  return(membership_trajectories)
}


theta_trajectories_summary <- function(VI_result){
  num_iter <- length(VI_result$theta)
  num_parameters <- 3
  theta_trajectories <- data.frame(Iter = rep(1:num_iter, each=num_parameters),
                                        Block = rep(c("Diag1", "OffDiag", "Diag2"), num_iter),
                                        Probability = 0)
  theta_trajectories$Probability <- do.call(c, lapply(VI_result$theta,
                                                      function(element) element[c(1,3,4)]))
  return(theta_trajectories)
}


membership_trajectories_degree <- membership_trajectories_summary(result_degree)
theta_trajectories_degree <- theta_trajectories_summary(result_degree)
membership_trajectories_sc <- membership_trajectories_summary(result_sc)
theta_trajectories_sc <- theta_trajectories_summary(result_sc)

membership_plot_degree <- ggplot(membership_trajectories_degree, aes(x=Iter, y=Membership, group=Taxon)) +
  geom_point(size=1.2, alpha=0.6) +
  geom_line(size=0.5, alpha=0.1) +
  theme(legend.position="none")

theta_plot_degree <- ggplot(theta_trajectories_degree, aes(x=Iter, y=Probability, group=Block))+
  geom_point(size=1.2, aes(color=Block), alpha=0.8) +
  geom_line(size=0.5, aes(color=Block), alpha=0.4)

membership_plot_sc <- ggplot(membership_trajectories_sc, aes(x=Iter, y=Membership, group=Taxon)) +
  geom_point(size=1.2, alpha=0.6) +
  geom_line(size=0.5, alpha=0.1) +
  theme(legend.position="none")

theta_plot_sc <- ggplot(theta_trajectories_sc, aes(x=Iter, y=Probability, group=Block))+
  geom_point(size=1.2, aes(color=Block), alpha=0.8) +
  geom_line(size=0.5, aes(color=Block), alpha=0.4)


library(cowplot)

plot_grid(membership_plot_degree, theta_plot_degree, membership_plot_sc, theta_plot_sc,
          nrow=2, rel_widths = c(1, 1.4), labels="AUTO")


# simplest SBM model with binary entries
taxaSBM <- BM_bernoulli(membership_type="SBM_sym", adj=GLM_decision_mat_high,
                        explore_min=2, explore_max=2, plotting='', verbosity=6)
taxaSBM$estimate()

result_pred_prob_package <- taxaSBM$memberships[[2]]$Z
difference <- 1 - result_pred_prob_sc - result_pred_prob_package
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
