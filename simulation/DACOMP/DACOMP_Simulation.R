rm(list=ls())
library(dplyr)
library(dacomp)

data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
dacomp_folder <- '/home/wangmk/MDAWG/POLDA/simulation/DACOMP'


source(file.path(dacomp_folder, "DACOMP_utils.R"))
ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 10


# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
dacomp_null <- pipeline(AGP_null)
performance_null <- evaluation(taxa_info_null, dacomp_null, nullcase=TRUE)
output_null <- list(performance=performance_null,
                    dacomp_result=dacomp_null)
saveRDS(output_null, 
        file.path(dacomp_folder, "null", 
                  sprintf("summary_null_%d.rds", ID)))

# rare DA taxa, low proportion(5%)
AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
                                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
dacomp_unbalanced_rare_low <- pipeline(AGP_unbalanced_rare_low)
performance_unbalanced_rare_low <- evaluation(taxa_info_unbalanced_rare_low,
                                              dacomp_unbalanced_rare_low)
output_unbalanced_rare_low <- list(performance=performance_unbalanced_rare_low,
                                   dacomp_result=dacomp_unbalanced_rare_low)
saveRDS(output_unbalanced_rare_low, 
        file.path(dacomp_folder, "DA_rare", "low", 
                  sprintf("summary_unbalanced_rare_low_%d.rds", ID)))


# rare DA taxa, medium proportion(10%)
AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                                sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
dacomp_unbalanced_rare_medium <- pipeline(AGP_unbalanced_rare_medium)
performance_unbalanced_rare_medium <- evaluation(taxa_info_unbalanced_rare_medium,
                                                 dacomp_unbalanced_rare_medium)
output_unbalanced_rare_medium <- list(performance=performance_unbalanced_rare_medium,
                                   dacomp_result=dacomp_unbalanced_rare_medium)
saveRDS(output_unbalanced_rare_medium, 
        file.path(dacomp_folder, "DA_rare", "medium", 
                  sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))


# rare DA taxa, high proportion(20%)
AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                                sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
dacomp_unbalanced_rare_high <- pipeline(AGP_unbalanced_rare_high)
performance_unbalanced_rare_high <- evaluation(taxa_info_unbalanced_rare_high,
                                                 dacomp_unbalanced_rare_high)
output_unbalanced_rare_high <- list(performance=performance_unbalanced_rare_high,
                                      dacomp_result=dacomp_unbalanced_rare_high)
saveRDS(output_unbalanced_rare_high, 
        file.path(dacomp_folder, "DA_rare", "high", 
                  sprintf("summary_unbalanced_rare_high_%d.rds", ID)))


# abundant DA taxa, low proportion(5%)
AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
                                                 sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
dacomp_unbalanced_abundant_low <- pipeline(AGP_unbalanced_abundant_low)
performance_unbalanced_abundant_low <- evaluation(taxa_info_unbalanced_abundant_low,
                                                  dacomp_unbalanced_abundant_low)
output_unbalanced_abundant_low <- list(performance=performance_unbalanced_abundant_low,
                                       dacomp_result=dacomp_unbalanced_abundant_low)
saveRDS(output_unbalanced_abundant_low, 
        file.path(dacomp_folder, "DA_abundant", "low", 
                  sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))


# abundant DA taxa, medium proportion(10%)
AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
                                                    sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
dacomp_unbalanced_abundant_medium <- pipeline(AGP_unbalanced_abundant_medium)
performance_unbalanced_abundant_medium <- evaluation(taxa_info_unbalanced_abundant_medium,
                                                     dacomp_unbalanced_abundant_medium)
output_unbalanced_abundant_medium <- list(performance=performance_unbalanced_abundant_medium,
                                          dacomp_result=dacomp_unbalanced_abundant_medium)
saveRDS(output_unbalanced_abundant_medium, 
        file.path(dacomp_folder, "DA_abundant", "medium", 
                  sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))


# abundant DA taxa, high proportion(20%)
AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
                                                  sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
dacomp_unbalanced_abundant_high <- pipeline(AGP_unbalanced_abundant_high)
performance_unbalanced_abundant_high <- evaluation(taxa_info_unbalanced_abundant_high,
                                                   dacomp_unbalanced_abundant_high)
output_unbalanced_abundant_high <- list(performance=performance_unbalanced_abundant_high,
                                        dacomp_result=dacomp_unbalanced_abundant_high)
saveRDS(output_unbalanced_abundant_high, 
        file.path(dacomp_folder, "DA_abundant", "high", 
                  sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))


