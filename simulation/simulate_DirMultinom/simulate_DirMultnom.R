
remove(list=ls())

source("simulation/simulate_DirMultinom/utils.R")

AGP_template <- readRDS('simulation/data/AGP_template.rds')
AGP_count_mat <- otu_table(AGP_template)@.Data

AGP_parameters <- estimate_param(AGP_count_mat)
AGP_dirgamma <- AGP_parameters$dirmult_composition

# the estimated gamma generates too much overdispersion, so I multiple all the parameters by the same constant
scaled_AGP_dirgamma <- AGP_dirgamma * 1 / min(AGP_dirgamma) # this need to change


# null sample
AGP_null <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                           seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                           seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                           nSample = 100,
                                                           grp_ratio=1,
                                                           DA_avg_logfold = 3,
                                                           DA_sd_logfold=0.1,
                                                           DA_proportion=0, seed=1, zinf=0.5,
                                                           DA_direction="balanced",
                                                           DA_mode="mixed")

# DA taxa include both rare and abundant ones, changes include both upward and downward
AGP_simulation_balanced_mixed <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                  seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                  seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                  nSample = 100,
                                  grp_ratio=1,
                                  DA_avg_logfold = 3,
                                  DA_sd_logfold=0.1,
                                  DA_proportion=0.2, seed=1, zinf=0.5,
                                  DA_direction="balanced",
                                  DA_mode="mixed")

taxa_info_balanced_mixed <- AGP_simulation_balanced_mixed$taxa_info

# DA taxa include only abundant ones, changes include both upward and downward
AGP_simulation_balanced_abundant <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                  seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                  seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                  nSample = 100,
                                                  grp_ratio=1,
                                                  DA_avg_logfold = 3,
                                                  DA_sd_logfold=0.1,
                                                  DA_proportion=0.2, seed=1,
                                                  DA_direction="balanced",
                                                  DA_mode="abundant")

taxa_info_balanced_abundant <- AGP_simulation_balanced_abundant$taxa_info

# DA taxa include both rare and abundant taxa, changes upward only
AGP_simulation_unbalanced_mixed <- SimulateCount(dirmult_composition = AGP_parameters$dirmult_composition,
                                                  seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                  seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                  nSample = 100,
                                                  grp_ratio=1,
                                                  DA_avg_logfold = 3,
                                                  DA_sd_logfold=0.1,
                                                  DA_proportion=0.2, seed=1,
                                                  DA_direction="unbalanced",
                                                  DA_mode="mixed")

taxa_info_unbalanced_mixed <- AGP_simulation_unbalanced_mixed$taxa_info

# DA taxa include only abundant taxa, changes upward only
AGP_simulation_unbalanced_abundant <- SimulateCount(dirmult_composition = AGP_parameters$dirmult_composition,
                                                 seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                 seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                 nSample = 100,
                                                 grp_ratio=1,
                                                 DA_avg_logfold = 3,
                                                 DA_sd_logfold=0.1,
                                                 DA_proportion=0.2, seed=1,
                                                 DA_direction="unbalanced",
                                                 DA_mode="abundant")

taxa_info_unbalanced_abundant <- AGP_simulation_unbalanced_abundant$taxa_info

# DA taxa include only rare taxa, changes upward only
AGP_simulation_unbalanced_rare <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                               seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                               seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                               nSample = 100,
                                               grp_ratio=1,
                                               DA_avg_logfold = 2,
                                               DA_sd_logfold=0.1,
                                               DA_proportion=0.2, seed=1, zinf=0.5,
                                               DA_direction="unbalanced",
                                               DA_mode="rare")


source("POLDA/POLDA.R") # parallel function need to change for POLDA.R


polda_result_balanced_mixed <- polda(otu_table=AGP_simulation_balanced_mixed$count_mat, 
                                     metadata=AGP_simulation_balanced_mixed$sample_metadata,
                      covar="X")

polda_result_unbalanced_mixed <- polda(otu_table=AGP_simulation_unbalanced_mixed$count_mat, 
                                       metadata=AGP_simulation_unbalanced_mixed$sample_metadata,
                                       covar="X")

polda_result_unbalanced_abundant <- polda(otu_table=AGP_simulation_unbalanced_abundant$count_mat, 
                                          metadata=AGP_simulation_unbalanced_abundant$sample_metadata,
                                          covar="X")

polda_result_unbalanced_rare <- polda(otu_table=AGP_simulation_unbalanced_rare$count_mat, 
                                      metadata=AGP_simulation_unbalanced_rare$sample_metadata,
                                      covar="X")


