
remove(list=ls())

source("simulation/simulate_DirMultinom/utils.R")

AGP_template <- readRDS('simulation/data/AGP_template.rds')
AGP_count_mat <- otu_table(AGP_template)@.Data

AGP_parameters <- estimate_param(AGP_count_mat, seed=3)
AGP_dirgamma <- AGP_parameters$dirmult_composition

# the estimated gamma generates too much overdispersion, so I multiple all the parameters by the same constant
scaled_AGP_dirgamma <- AGP_dirgamma * 5 


# null sample
AGP_null <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                            seqdepth_mean = AGP_parameters$Seqdepth_mean,
                            seqdepth_theta = AGP_parameters$Seqdepth_theta,
                            nSample=100,
                            DA_proportion=0, seed=2, zinf=0.7)


# evaluate POLDA performance
evaluation <- function(taxa_truth, polda_result){
  true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
  check_reftaxa <- polda_result$Reference_Taxa %in% true_DA_taxa
  all_DA_taxa <- c(polda_result$Structural_zero_Taxa, polda_result$DA_taxa)
  check_DAtaxa <- all_DA_taxa %in% true_DA_taxa
  reftaxa_error <- mean(check_reftaxa)
  FDR <- 1 - mean(check_DAtaxa)
  Power <- sum(check_DAtaxa) / length(true_DA_taxa)
  result <- list(reftaxa_error=reftaxa_error,
                 FDR = FDR,
                 Power=Power)
  return(result)
}

source("POLDA/POLDA.R") # parallel function need to change for POLDA.R

polda_null <- polda(otu_table=AGP_null$count_mat,
                      metadata=AGP_null$sample_metadata,
                      covar="X")

# DA taxa include only rare taxa, changes upward only
AGP_simulation_unbalanced_rare <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                nSample = 100,
                                                grp_ratio=1,
                                                DA_avg_logfold = 2,
                                                DA_sd_logfold=0.1,
                                                DA_proportion=0.1, seed=28, zinf=0.6,
                                                DA_direction="unbalanced",
                                                DA_mode="rare")

polda_result_unbalanced_rare <- polda(otu_table=AGP_simulation_unbalanced_rare$count_mat, 
                                      metadata=AGP_simulation_unbalanced_rare$sample_metadata,
                                      covar="X")

taxa_info_unbalanced_rare <- AGP_simulation_unbalanced_rare$taxa_info

performance_unbalanced_rare <- evaluation(taxa_truth = taxa_info_unbalanced_rare,
                                          polda_result = polda_result_unbalanced_rare)

performance_unbalanced_rare

# DA taxa include only abundant taxa, changes upward only
AGP_simulation_unbalanced_abundant <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma,
                                                    seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                    seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                    nSample = 100,
                                                    grp_ratio=1,
                                                    DA_avg_logfold = 1,
                                                    DA_sd_logfold=0.1,
                                                    DA_proportion=0.2, seed=26, zinf=0.6,
                                                    DA_direction="unbalanced",
                                                    DA_mode="abundant")

taxa_info_unbalanced_abundant <- AGP_simulation_unbalanced_abundant$taxa_info

polda_result_unbalanced_abundant <- polda(otu_table=AGP_simulation_unbalanced_abundant$count_mat, 
                                          metadata=AGP_simulation_unbalanced_abundant$sample_metadata,
                                          covar="X")

performance_unbalanced_abundant <- evaluation(taxa_truth = taxa_info_unbalanced_abundant,
                                              polda_result = polda_result_unbalanced_abundant)

performance_unbalanced_abundant



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


