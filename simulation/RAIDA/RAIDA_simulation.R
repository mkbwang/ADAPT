rm(list=ls())
library(protoclust)
library(limma)
source('/home/wangmk/MDAWG/RAIDA/RAIDA/R/raida.R')


data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
RAIDA_folder <- '/home/wangmk/MDAWG/POLDA/simulation/RAIDA'
source(file.path(RAIDA_folder, 'RAIDA_utils.R'))

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 10

n.samples <- c(50, 50)
# null
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata_null <- AGP_null$sample_metadata
count_mat_null <- reorganize(AGP_null)
raida_null <- raida(count_mat_null, n.samples, 
                    show.ref.features = TRUE, show.all.features = TRUE)
performance_null <- evaluation(taxa_info_null, raida_null, nullcase=TRUE)
output_null <- list(performance=performance_null,
                    raida_result=raida_null)
saveRDS(output_null, 
        file.path(RAIDA_folder, "null", 
                  sprintf("summary_null_%d.rds", ID)))



# rare DA taxa, low proportion(5%)
AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
                                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
sample_metadata_rare_low <- AGP_unbalanced_rare_low$sample_metadata
count_mat_rare_low <- reorganize(AGP_unbalanced_rare_low)
raida_rare_low <- raida(count_mat_rare_low, n.samples, 
                        show.ref.features = TRUE, show.all.features = TRUE)
performance_rare_low <- evaluation(taxa_info_unbalanced_rare_low, raida_rare_low)
output_rare_low <- list(performance=performance_rare_low,
                    raida_result=raida_rare_low)
saveRDS(output_rare_low, 
        file.path(RAIDA_folder, "DA_rare", "low", 
                  sprintf("summary_unbalanced_rare_low_%d.rds", ID)))


# rare DA taxa, medium proportion(10%)
AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                             sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
sample_metadata_rare_medium <- AGP_unbalanced_rare_medium$sample_metadata
count_mat_rare_medium <- reorganize(AGP_unbalanced_rare_medium)
raida_rare_medium <- raida(count_mat_rare_medium, n.samples, 
                        show.ref.features = TRUE, show.all.features = TRUE)
performance_rare_medium <- evaluation(taxa_info_unbalanced_rare_medium, raida_rare_medium)
output_rare_medium <- list(performance=performance_rare_medium,
                        raida_result=raida_rare_medium)
saveRDS(output_rare_medium, 
        file.path(RAIDA_folder, "DA_rare", "medium", 
                  sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))


# rare DA taxa, high proportion(20%)
AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                                sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
sample_metadata_rare_high <- AGP_unbalanced_rare_high$sample_metadata
count_mat_rare_high <- reorganize(AGP_unbalanced_rare_high)
raida_rare_high <- raida(count_mat_rare_high, n.samples, 
                           show.ref.features = TRUE, show.all.features = TRUE)
performance_rare_high <- evaluation(taxa_info_unbalanced_rare_medium, raida_rare_high)
output_rare_high <- list(performance=performance_rare_high,
                           raida_result=raida_rare_high)
saveRDS(output_rare_high, 
        file.path(RAIDA_folder, "DA_rare", "high", 
                  sprintf("summary_unbalanced_rare_high_%d.rds", ID)))


# abundant DA taxa, low proportion(5%)
AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
                                                 sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
sample_metadata_abundant_low <- AGP_unbalanced_abundant_low$sample_metadata
count_mat_abundant_low <- reorganize(AGP_unbalanced_abundant_low)
raida_abundant_low <- raida(count_mat_abundant_low, n.samples, 
                            show.ref.features = TRUE, show.all.features = TRUE)
performance_abundant_low <- evaluation(taxa_info_unbalanced_abundant_low, raida_abundant_low)
output_abundant_low <- list(performance=performance_abundant_low,
                         raida_result=raida_abundant_low)
saveRDS(output_abundant_low, 
        file.path(RAIDA_folder, "DA_abundant", "low", 
                  sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))


# abundant DA taxa, medium proportion(10%)
AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
                                                    sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
sample_metadata_abundant_medium <- AGP_unbalanced_abundant_medium$sample_metadata
count_mat_abundant_medium <- reorganize(AGP_unbalanced_abundant_medium)
raida_abundant_medium <- raida(count_mat_abundant_medium, n.samples, 
                               show.ref.features = TRUE, show.all.features = TRUE)
performance_abundant_medium <- evaluation(taxa_info_unbalanced_abundant_medium, raida_abundant_medium)
output_abundant_medium <- list(performance=performance_abundant_medium,
                            raida_result=raida_abundant_medium)
saveRDS(output_abundant_medium, 
        file.path(RAIDA_folder, "DA_abundant", "medium", 
                  sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))


# abundant DA taxa, high proportion(20%)
AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
                                                  sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
sample_metadata_abundant_high <- AGP_unbalanced_abundant_high$sample_metadata
count_mat_abundant_high <- reorganize(AGP_unbalanced_abundant_high)
raida_abundant_high <- raida(count_mat_abundant_high, n.samples, 
                             show.ref.features = TRUE, show.all.features = TRUE)
performance_abundant_high <- evaluation(taxa_info_unbalanced_abundant_medium, raida_abundant_high)
output_abundant_high <- list(performance=performance_abundant_high,
                               raida_result=raida_abundant_high)
saveRDS(output_abundant_high, 
        file.path(RAIDA_folder, "DA_abundant", "high", 
                  sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))






