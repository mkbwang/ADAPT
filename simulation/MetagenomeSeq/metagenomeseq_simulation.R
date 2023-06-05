rm(list=ls())
library(dplyr)
library(metagenomeSeq)

data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
metagenomeseq_folder <- '/home/wangmk/MDAWG/POLDA/simulation/MetagenomeSeq'

source(file.path(metagenomeseq_folder, "metagenomeseq_utils.R"))
ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 10

# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata_null <- AGP_null$sample_metadata
count_mat_null <- AGP_null$count_mat
MRobj_null <- setExperiment(sample_metadata_null, count_mat_null)
# MRobj_null <- wrenchNorm(MRobj_null, condition=MRobj_null$X)
mod <- model.matrix(~X, data = sample_metadata_null)
MRNull_res <- fitFeatureModel(MRobj_null, mod)
summary_null <- evaluation(taxa_info_null, MRNull_res, nullcase=TRUE)
saveRDS(summary_null, 
        file.path(metagenomeseq_folder, "null", 
                  sprintf("summary_null_%d.rds", ID)))



# rare DA taxa, low proportion(5%)
AGP_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
                                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
taxa_info_rare_low <- AGP_rare_low$taxa_info
sample_metadata_rare_low <- AGP_rare_low$sample_metadata
count_mat_rare_low <- AGP_rare_low$count_mat
MRobj_rare_low <- setExperiment(sample_metadata_rare_low, count_mat_rare_low)
mod <- model.matrix(~X, data = sample_metadata_rare_low)
MR_rare_low_res <- fitFeatureModel(MRobj_rare_low, mod)
summary_rare_low <- evaluation(taxa_info_rare_low,
                                             MR_rare_low_res)
saveRDS(summary_rare_low, 
        file.path(metagenomeseq_folder, "DA_rare", "low", 
                  sprintf("summary_unbalanced_rare_low_%d.rds", ID)))


# rare DA taxa, medium proportion(10%)
AGP_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                  sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
taxa_info_rare_medium <- AGP_rare_medium$taxa_info
sample_metadata_rare_medium <- AGP_rare_medium$sample_metadata
count_mat_rare_medium <- AGP_rare_medium$count_mat
MRobj_rare_medium <- setExperiment(sample_metadata_rare_medium, count_mat_rare_medium)
mod <- model.matrix(~X, data = sample_metadata_rare_medium)
MR_rare_medium_res <- fitFeatureModel(MRobj_rare_medium, mod)
summary_rare_medium <- evaluation(taxa_info_rare_medium,
                               MR_rare_medium_res)
saveRDS(summary_rare_medium, 
        file.path(metagenomeseq_folder, "DA_rare", "medium", 
                  sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))



# rare DA taxa, high proportion(20%)
AGP_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                     sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
taxa_info_rare_high <- AGP_rare_high$taxa_info
sample_metadata_rare_high <- AGP_rare_high$sample_metadata
count_mat_rare_high <- AGP_rare_high$count_mat
MRobj_rare_high <- setExperiment(sample_metadata_rare_high, count_mat_rare_high)
mod <- model.matrix(~X, data = sample_metadata_rare_high)
MR_rare_high_res <- fitFeatureModel(MRobj_rare_high, mod)
summary_rare_high <- evaluation(taxa_info_rare_high,
                                  MR_rare_high_res)
saveRDS(summary_rare_high, 
        file.path(metagenomeseq_folder, "DA_rare", "high", 
                  sprintf("summary_unbalanced_rare_high_%d.rds", ID)))



# abundant DA taxa, low proportion(5%)
AGP_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
                                      sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
taxa_info_abundant_low <- AGP_abundant_low$taxa_info
sample_metadata_abundant_low <- AGP_abundant_low$sample_metadata
count_mat_abundant_low <- AGP_abundant_low$count_mat
MRobj_abundant_low <- setExperiment(sample_metadata_abundant_low, count_mat_abundant_low)
mod <- model.matrix(~X, data = sample_metadata_abundant_low)
MR_abundant_low_res <- fitFeatureModel(MRobj_abundant_low, mod)
summary_abundant_low <- evaluation(taxa_info_abundant_low,
                                   MR_abundant_low_res)
saveRDS(summary_abundant_low, 
        file.path(metagenomeseq_folder, "DA_abundant", "low", 
                  sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))


# abundant DA taxa, medium proportion(10%)
AGP_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
                                         sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
taxa_info_abundant_medium <- AGP_abundant_medium$taxa_info
sample_metadata_abundant_medium <- AGP_abundant_medium$sample_metadata
count_mat_abundant_medium <- AGP_abundant_medium$count_mat
MRobj_abundant_medium <- setExperiment(sample_metadata_abundant_medium, count_mat_abundant_medium)
mod <- model.matrix(~X, data = sample_metadata_abundant_medium)
MR_abundant_medium_res <- fitFeatureModel(MRobj_abundant_medium, mod)
summary_abundant_medium <- evaluation(taxa_info_abundant_medium,
                                      MR_abundant_medium_res)
saveRDS(summary_abundant_medium, 
        file.path(metagenomeseq_folder, "DA_abundant", "medium", 
                  sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))



# abundant DA taxa, high proportion(20%)
AGP_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
                                       sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
taxa_info_abundant_high <- AGP_abundant_high$taxa_info
sample_metadata_abundant_high <- AGP_abundant_high$sample_metadata
count_mat_abundant_high <- AGP_abundant_high$count_mat
MRobj_abundant_high <- setExperiment(sample_metadata_abundant_high, count_mat_abundant_high)
mod <- model.matrix(~X, data = sample_metadata_abundant_high)
MR_abundant_high_res <- fitFeatureModel(MRobj_abundant_high, mod)
summary_abundant_high <- evaluation(taxa_info_abundant_high,
                                    MR_abundant_high_res)
saveRDS(summary_abundant_high, 
        file.path(metagenomeseq_folder, "DA_abundant", "high", 
                  sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))



