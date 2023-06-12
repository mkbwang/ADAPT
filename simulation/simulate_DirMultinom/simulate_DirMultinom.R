library(ggplot2)
remove(list=ls())

folder <- '/home/wangmk/MDAWG/POLDA/simulation'
source(file.path(folder, "simulate_DirMultinom/utils.R"))

# AGP_template <- readRDS('simulation/data/AGP_template.rds')
# AGP_count_mat <- otu_table(AGP_template)@.Data
# 
# AGP_parameters <- estimate_param(AGP_count_mat, seed=5)
# AGP_dirgamma <- AGP_parameters$dirmult_composition
# scaled_AGP_dirgamma <- AGP_dirgamma * 5
# mean_depth <- AGP_parameters$Seqdepth_mean
# theta_depth <- AGP_parameters$Seqdepth_theta
# save(scaled_AGP_dirgamma, mean_depth, theta_depth, file="simulation/data/AGP_dirgamma.RData")
# the estimated gamma generates too much overdispersion, so I multiple all the parameters by the same constant
load(file.path(folder, "data/AGP_dirgamma.RData"))
library(cowplot)
example_composition <- rdirichlet(scaled_AGP_dirgamma, nSample=2, seed=10)
expected_composition <- scaled_AGP_dirgamma / sum(scaled_AGP_dirgamma)
expected_composition_df <- data.frame(composition=expected_composition)
expected_composition_plot <- ggplot(expected_composition_df, aes(x=composition))+
        geom_histogram(bins=20, color="black", fill="white")+
        scale_x_log10()+
        xlab("") +ylab("Frequency") + 
        ggtitle("Expected Compositions of 1000 Taxa")+theme_bw()+
        theme(text=element_text(size=12),
              plot.title = element_text(hjust = 0.5))

sample_depths <- rnegbin(1000, mu=mean_depth, theta=theta_depth)
sample_depths_df <- data.frame(Depths=sample_depths)
sample_depths_plot <- ggplot(sample_depths_df, aes(x=Depths)) + 
        geom_histogram(binwidth=10000, color="black", fill="white")+
        xlab("") +ylab("Frequency") + 
        ggtitle("Distribution of Simulated Sequencing Depths")+theme_bw()+
        theme(text=element_text(size=12),
              plot.title = element_text(hjust = 0.5))

plot_grid(expected_composition_plot,
          sample_depths_plot, ncol=2)

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# null sample
# AGP_null <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
#                             seqdepth_mean = AGP_parameters$Seqdepth_mean,
#                             seqdepth_theta = AGP_parameters$Seqdepth_theta,
#                             nSample=100,
#                             DA_proportion=0, seed=2, zinf=0.6)



# evaluate POLDA performance
# evaluation <- function(taxa_truth, polda_result){
#   true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
#   check_reftaxa <- polda_result$Reference_Taxa %in% true_DA_taxa
#   all_DA_taxa <- c(polda_result$Structural_zero_Taxa, polda_result$DA_taxa)
#   check_DAtaxa <- all_DA_taxa %in% true_DA_taxa
#   reftaxa_error <- mean(check_reftaxa)
#   FDR <- 1 - mean(check_DAtaxa)
#   Power <- sum(check_DAtaxa) / length(true_DA_taxa)
#   result <- list(reftaxa_error=reftaxa_error,
#                  FDR = FDR,
#                  Power=Power)
#   return(result)
# }

# source("POLDA/POLDA.R") # parallel function need to change for POLDA.R

# polda_null <- polda(otu_table=AGP_null$count_mat,
#                       metadata=AGP_null$sample_metadata,
#                       covar="X")

# DA taxa include only rare taxa, changes upward only
# low proportion: 5%
AGP_simulation_unbalanced_rare_low <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                seqdepth_mean = mean_depth,
                                                seqdepth_theta = theta_depth,
                                                nSample = 100,
                                                grp_ratio=1,
                                                DA_avg_logfold = 2,
                                                DA_sd_logfold=0.1,
                                                DA_proportion=0.05, seed=ID, zinf=0.6,
                                                DA_direction="unbalanced",
                                                DA_mode="rare")
saveRDS(AGP_simulation_unbalanced_rare_low, 
        file.path(folder, "data", "DA_rare", "low", 
                  sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))

example <- readRDS(file.path(folder, "data", "DA_rare", "low", 
                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", 20)))

# medium proportion: 10%
AGP_simulation_unbalanced_rare_medium <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                    seqdepth_mean = mean_depth,
                                                    seqdepth_theta = theta_depth,
                                                    nSample = 100,
                                                    grp_ratio=1,
                                                    DA_avg_logfold = 2,
                                                    DA_sd_logfold=0.1,
                                                    DA_proportion=0.1, seed=ID, zinf=0.6,
                                                    DA_direction="unbalanced",
                                                    DA_mode="rare")
saveRDS(AGP_simulation_unbalanced_rare_medium, 
        file.path(folder, "data", "DA_rare", "medium", 
                  sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))

example <- readRDS(file.path(folder, "data", "DA_rare", "medium", 
                             sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", 20)))

# high proportion: 20%
AGP_simulation_unbalanced_rare_high <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                       seqdepth_mean = mean_depth,
                                                       seqdepth_theta = theta_depth,
                                                       nSample = 100,
                                                       grp_ratio=1,
                                                       DA_avg_logfold = 2,
                                                       DA_sd_logfold=0.1,
                                                       DA_proportion=0.2, seed=ID, zinf=0.6,
                                                       DA_direction="unbalanced",
                                                       DA_mode="rare")
saveRDS(AGP_simulation_unbalanced_rare_high, 
        file.path(folder, "data", "DA_rare", "high", 
                  sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))

example <- readRDS(file.path(folder, "data", "DA_rare", "high", 
                             sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", 20)))


# DA taxa include only abundant taxa, changes upward only
# low proportion: 5%
AGP_simulation_unbalanced_abundant_low <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                    seqdepth_mean = mean_depth,
                                                    seqdepth_theta = theta_depth,
                                                    nSample = 100,
                                                    grp_ratio=1,
                                                    DA_avg_logfold = 1,
                                                    DA_sd_logfold=0.1,
                                                    DA_proportion=0.05, seed=ID, zinf=0.6,
                                                    DA_direction="unbalanced",
                                                    DA_mode="abundant")
saveRDS(AGP_simulation_unbalanced_abundant_low, 
        file.path(folder, "data", "DA_abundant", "low", 
                  sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))


# medium proportion: 10%
AGP_simulation_unbalanced_abundant_medium <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                        seqdepth_mean = mean_depth,
                                                        seqdepth_theta = theta_depth,
                                                        nSample = 100,
                                                        grp_ratio=1,
                                                        DA_avg_logfold = 1,
                                                        DA_sd_logfold=0.1,
                                                        DA_proportion=0.1, seed=ID, zinf=0.6,
                                                        DA_direction="unbalanced",
                                                        DA_mode="abundant")
saveRDS(AGP_simulation_unbalanced_abundant_medium, 
        file.path(folder, "data", "DA_abundant", "medium", 
                  sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))


# high proportion: 20%
AGP_simulation_unbalanced_abundant_high <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                           seqdepth_mean = mean_depth,
                                                           seqdepth_theta = theta_depth,
                                                           nSample = 100,
                                                           grp_ratio=1,
                                                           DA_avg_logfold = 1,
                                                           DA_sd_logfold=0.1,
                                                           DA_proportion=0.2, seed=ID, zinf=0.6,
                                                           DA_direction="unbalanced",
                                                           DA_mode="abundant")
saveRDS(AGP_simulation_unbalanced_abundant_high, 
        file.path(folder, "data", "DA_abundant", "high", 
                  sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))


# # DA taxa include both rare and abundant ones, changes include both upward and downward
# AGP_simulation_balanced_mixed <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
#                                   seqdepth_mean = AGP_parameters$Seqdepth_mean,
#                                   seqdepth_theta = AGP_parameters$Seqdepth_theta,
#                                   nSample = 100,
#                                   grp_ratio=1,
#                                   DA_avg_logfold = 2,
#                                   DA_sd_logfold=0.1,
#                                   DA_proportion=0.10, seed=7, zinf=0.7,
#                                   DA_direction="balanced",
#                                   DA_mode="mixed")
# 
# taxa_info_balanced_mixed <- AGP_simulation_balanced_mixed$taxa_info
# polda_result_balanced_mixed <- polda(otu_table=AGP_simulation_balanced_mixed$count_mat, 
#                                      metadata=AGP_simulation_balanced_mixed$sample_metadata,
#                                      covar="X")
# 
# performance_balanced_mixed <- evaluation(taxa_truth = taxa_info_balanced_mixed,
#                                          polda_result = polda_result_balanced_mixed)
# performance_balanced_mixed
# 
# # DA taxa include only abundant ones, changes include both upward and downward
# AGP_simulation_balanced_abundant <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
#                                                   seqdepth_mean = AGP_parameters$Seqdepth_mean,
#                                                   seqdepth_theta = AGP_parameters$Seqdepth_theta,
#                                                   nSample = 100,
#                                                   grp_ratio=1,
#                                                   DA_avg_logfold = 2,
#                                                   DA_sd_logfold=0.1,
#                                                   DA_proportion=0.1, seed=80, zinf=0.7,
#                                                   DA_direction="balanced",
#                                                   DA_mode="abundant")
# 
# taxa_info_balanced_abundant <- AGP_simulation_balanced_abundant$taxa_info
# polda_result_balanced_abundant <- polda(otu_table=AGP_simulation_balanced_abundant$count_mat, 
#                                         metadata=AGP_simulation_balanced_abundant$sample_metadata,
#                                         covar="X")
# 
# 
# performance_balanced_abundant <- evaluation(taxa_truth = taxa_info_balanced_abundant,
#                                          polda_result = polda_result_balanced_abundant)
# 
# performance_balanced_abundant
# 
# # DA taxa include only rare taxa, changes both upward and downward
# AGP_simulation_balanced_rare <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
#                                                   seqdepth_mean = AGP_parameters$Seqdepth_mean,
#                                                   seqdepth_theta = AGP_parameters$Seqdepth_theta,
#                                                   nSample = 100,
#                                                   grp_ratio=1,
#                                                   DA_avg_logfold = 2,
#                                                   DA_sd_logfold=0.1,
#                                                   DA_proportion=0.1, seed=6, zinf=0.7,
#                                                   DA_direction="balanced",
#                                                   DA_mode="rare")
# 
# taxa_info_balanced_rare <- AGP_simulation_balanced_rare$taxa_info
# polda_result_balanced_rare <- polda(otu_table=AGP_simulation_balanced_rare$count_mat, 
#                                         metadata=AGP_simulation_balanced_rare$sample_metadata,
#                                         covar="X")
# 
# performance_balanced_rare <- evaluation(taxa_truth = taxa_info_balanced_rare,
#                                          polda_result = polda_result_balanced_rare)
# 
# performance_balanced_rare
# 
# 
# # DA taxa include both abundant and rare taxa, changes upward only
# AGP_simulation_unbalanced_mixed <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
#                                                  seqdepth_mean = AGP_parameters$Seqdepth_mean,
#                                                  seqdepth_theta = AGP_parameters$Seqdepth_theta,
#                                                  nSample = 100,
#                                                  grp_ratio=1,
#                                                  DA_avg_logfold = 2,
#                                                  DA_sd_logfold=0.1,
#                                                  DA_proportion=0.2, seed=20, zinf=0.7,
#                                                  DA_direction="unbalanced",
#                                                  DA_mode="mixed")
# 
# taxa_info_unbalanced_mixed <- AGP_simulation_unbalanced_mixed$taxa_info
# polda_result_unbalanced_mixed <- polda(otu_table=AGP_simulation_unbalanced_mixed$count_mat, 
#                                     metadata=AGP_simulation_unbalanced_mixed$sample_metadata,
#                                     covar="X")
# 
# performance_unbalanced_mixed <- evaluation(taxa_truth = taxa_info_unbalanced_mixed,
#                                         polda_result = polda_result_unbalanced_mixed)
# 
# performance_unbalanced_mixed
# 
# 


