
remove(list=ls())

source("simulation/simulate_DirMultinom/utils.R")

AGP_template <- readRDS('simulation/data/AGP_template.rds')
AGP_count_mat <- otu_table(AGP_template)@.Data

AGP_parameters <- estimate_param(AGP_count_mat, seed=3)
AGP_dirgamma <- AGP_parameters$dirmult_composition

# the estimated gamma generates too much overdispersion, so I multiple all the parameters by the same constant
scaled_AGP_dirgamma1 <- AGP_dirgamma * 10 # this need to change
scaled_AGP_dirgamma2 <- AGP_dirgamma * 20


# null sample
AGP_null_0 <- SimulateCount(dirmult_composition = AGP_dirgamma,
                            seqdepth_mean = AGP_parameters$Seqdepth_mean,
                            seqdepth_theta = AGP_parameters$Seqdepth_theta,
                            DA_proportion=0, seed=1, zinf=0.5,
                            DA_direction="balanced",
                            DA_mode="mixed")

AGP_null_1 <- SimulateCount(dirmult_composition = AGP_dirgamma*5,
                            seqdepth_mean = AGP_parameters$Seqdepth_mean,
                            seqdepth_theta = AGP_parameters$Seqdepth_theta,
                            DA_proportion=0, seed=1, zinf=0.5,
                            DA_direction="balanced",
                            DA_mode="mixed")


AGP_null_2 <- SimulateCount(dirmult_composition = scaled_AGP_dirgamma2,
                            seqdepth_mean = AGP_parameters$Seqdepth_mean,
                            seqdepth_theta = AGP_parameters$Seqdepth_theta,
                            DA_proportion=0, seed=1, zinf=0.5,
                            DA_direction="balanced",
                            DA_mode="mixed")


evaluation <- function(taxa_truth, polda_result){
  true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
  check_reftaxa <- polda_result$Reference_Taxa %in% true_DA_taxa
  all_DA_taxa <- c(polda_result$DA_taxa, polda_result$Structural_zero_Taxa)
  check_DAtaxa <- all_DA_taxa %in% true_DA_taxa
  reftaxa_error <- mean(check_reftaxa)
  FDR <- 1-mean(check_DAtaxa)
  Power <- sum(check_DAtaxa) / length(true_DA_taxa)
  result <- list(reftaxa_error=reftaxa_error,
                 FDR = FDR,
                 Power=Power)
  return(result)
}

source("POLDA/POLDA.R") # parallel function need to change for POLDA.R

polda_null_0 <- polda(otu_table=AGP_null_0$count_mat,
                      metadata=AGP_null_0$sample_metadata,
                      covar="X")

hist(polda_null_0$P_Value$effect, nclass=40)

polda_null_1 <- polda(otu_table=AGP_null_1$count_mat,
                      metadata=AGP_null_1$sample_metadata,
                      covar="X")

hist(polda_null_1$P_Value$effect, nclass=40)

polda_null_2 <- polda(otu_table=AGP_null_2$count_mat,
                      metadata=AGP_null_2$sample_metadata,
                      covar="X")

hist(polda_null_2$P_Value$effect, nclass=40)

# DA taxa include both rare and abundant ones, changes include both upward and downward
AGP_simulation_balanced_mixed <- SimulateCount(dirmult_composition = AGP_dirgamma*5,
                                  seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                  seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                  nSample = 100,
                                  grp_ratio=1,
                                  DA_avg_logfold = 2,
                                  DA_sd_logfold=0.1,
                                  DA_proportion=0.2, seed=5, zinf=0.7,
                                  DA_direction="balanced",
                                  DA_mode="mixed")

taxa_info_balanced_mixed <- AGP_simulation_balanced_mixed$taxa_info
polda_result_balanced_mixed <- polda(otu_table=AGP_simulation_balanced_mixed$count_mat, 
                                     metadata=AGP_simulation_balanced_mixed$sample_metadata,
                                     covar="X")

performance_balanced_mixed <- evaluation(taxa_truth = taxa_info_balanced_mixed,
                                         polda_result = polda_result_balanced_mixed)
performance_balanced_mixed

# DA taxa include only abundant ones, changes include both upward and downward
AGP_simulation_balanced_abundant <- SimulateCount(dirmult_composition = AGP_dirgamma*5,
                                                  seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                  seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                  nSample = 100,
                                                  grp_ratio=1,
                                                  DA_avg_logfold = 2,
                                                  DA_sd_logfold=0.1,
                                                  DA_proportion=0.2, seed=5, zinf=0.7,
                                                  DA_direction="balanced",
                                                  DA_mode="abundant")

taxa_info_balanced_abundant <- AGP_simulation_balanced_abundant$taxa_info
polda_result_balanced_abundant <- polda(otu_table=AGP_simulation_balanced_abundant$count_mat, 
                                        metadata=AGP_simulation_balanced_abundant$sample_metadata,
                                        covar="X")


performance_balanced_abundant <- evaluation(taxa_truth = taxa_info_balanced_abundant,
                                         polda_result = polda_result_balanced_abundant)

performance_balanced_abundant

# DA taxa include only rare taxa, changes both upward and downward
AGP_simulation_balanced_rare <- SimulateCount(dirmult_composition = AGP_dirgamma*5,
                                                  seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                  seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                  nSample = 100,
                                                  grp_ratio=1,
                                                  DA_avg_logfold = 2,
                                                  DA_sd_logfold=0.1,
                                                  DA_proportion=0.2, seed=4, zinf=0.7,
                                                  DA_direction="balanced",
                                                  DA_mode="rare")

taxa_info_balanced_rare <- AGP_simulation_balanced_rare$taxa_info
polda_result_balanced_rare <- polda(otu_table=AGP_simulation_balanced_rare$count_mat, 
                                        metadata=AGP_simulation_balanced_rare$sample_metadata,
                                        covar="X")

performance_balanced_rare <- evaluation(taxa_truth = taxa_info_balanced_rare,
                                         polda_result = polda_result_balanced_rare)

performance_balanced_rare


# DA taxa include both abundant and rare taxa, changes upward only
AGP_simulation_unbalanced_mixed <- SimulateCount(dirmult_composition = AGP_dirgamma*5,
                                                 seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                 seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                 nSample = 100,
                                                 grp_ratio=1,
                                                 DA_avg_logfold = 2,
                                                 DA_sd_logfold=0.1,
                                                 DA_proportion=0.2, seed=20, zinf=0.7,
                                                 DA_direction="unbalanced",
                                                 DA_mode="mixed")

taxa_info_unbalanced_mixed <- AGP_simulation_unbalanced_mixed$taxa_info
polda_result_unbalanced_mixed <- polda(otu_table=AGP_simulation_unbalanced_mixed$count_mat, 
                                    metadata=AGP_simulation_unbalanced_mixed$sample_metadata,
                                    covar="X")

performance_unbalanced_mixed <- evaluation(taxa_truth = taxa_info_unbalanced_mixed,
                                        polda_result = polda_result_unbalanced_mixed)

performance_unbalanced_mixed

# DA taxa include only rare taxa, changes upward only
AGP_simulation_unbalanced_rare <- SimulateCount(dirmult_composition = AGP_dirgamma*5,
                                               seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                               seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                               nSample = 100,
                                               grp_ratio=1,
                                               DA_avg_logfold = 2,
                                               DA_sd_logfold=0.1,
                                               DA_proportion=0.2, seed=1, zinf=0.7,
                                               DA_direction="unbalanced",
                                               DA_mode="rare")

polda_result_unbalanced_rare <- polda(otu_table=AGP_simulation_unbalanced_rare$count_mat, 
                                       metadata=AGP_simulation_unbalanced_rare$sample_metadata,
                                       covar="X")

taxa_info_unbalanced_rare <- AGP_simulation_unbalanced_rare$taxa_info

performance_unbalanced_rare <- evaluation(taxa_truth = taxa_info_unbalanced_rare,
                                           polda_result = polda_result_unbalanced_rare)


# DA taxa include only abundant taxa, changes upward only
AGP_simulation_unbalanced_abundant <- SimulateCount(dirmult_composition = AGP_dirgamma*5,
                                                seqdepth_mean = AGP_parameters$Seqdepth_mean,
                                                seqdepth_theta = AGP_parameters$Seqdepth_theta,
                                                nSample = 100,
                                                grp_ratio=1,
                                                DA_avg_logfold = 2,
                                                DA_sd_logfold=0.1,
                                                DA_proportion=0.2, seed=1, zinf=0.7,
                                                DA_direction="unbalanced",
                                                DA_mode="abundant")

taxa_info_unbalanced_abundant <- AGP_simulation_unbalanced_abundant$taxa_info

polda_result_unbalanced_abundant <- polda(otu_table=AGP_simulation_unbalanced_abundant$count_mat, 
                                          metadata=AGP_simulation_unbalanced_abundant$sample_metadata,
                                          covar="X")

performance_unbalanced_abundant <- evaluation(taxa_truth = taxa_info_unbalanced_abundant,
                                          polda_result = polda_result_unbalanced_abundant)


