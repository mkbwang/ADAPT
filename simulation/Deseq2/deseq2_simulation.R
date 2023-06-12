rm(list=ls())
library(GUniFrac)
library(DESeq2)

data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
deseq2_folder <- '/home/wangmk/MDAWG/POLDA/simulation/Deseq2'

source(file.path(deseq2_folder, 'deseq2_utils.R'))

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 9

# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata_null <- AGP_null$sample_metadata
sample_metadata_null$X <- as.factor(sample_metadata_null$X)
count_mat_null <- AGP_null$count_mat
deseq_null <- gmpr_deseq(count_mat_null, sample_metadata_null)
summary_null <- evaluation(taxa_info_null, deseq_null, nullcase=T)
saveRDS(summary_null,
        file.path(deseq2_folder, "null",
                  sprintf("summary_null_%d.rds", ID)))


# rare DA taxa, low proportion(5%)
AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
                                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
sample_metadata_rare_low <- AGP_unbalanced_rare_low$sample_metadata
count_mat_rare_low <- AGP_unbalanced_rare_low$count_mat
deseq_rare_low <- gmpr_deseq(count_mat_rare_low, sample_metadata_rare_low)
summary_rare_low <- evaluation(taxa_info_unbalanced_rare_low, deseq_rare_low)
saveRDS(summary_rare_low,
        file.path(deseq2_folder, "DA_rare", "low",
                  sprintf("summary_unbalanced_rare_low_%d.rds", ID)))


# rare DA taxa, medium proportion(10%)
AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                                sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
sample_metadata_rare_medium <- AGP_unbalanced_rare_medium$sample_metadata
count_mat_rare_medium <- AGP_unbalanced_rare_medium$count_mat
deseq_rare_medium <- gmpr_deseq(count_mat_rare_medium, sample_metadata_rare_medium)
summary_rare_medium <- evaluation(taxa_info_unbalanced_rare_medium, deseq_rare_medium)
saveRDS(summary_rare_medium,
        file.path(deseq2_folder, "DA_rare", "medium",
                  sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))


# rare DA taxa, high proportion(20%)
AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                              sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
sample_metadata_rare_high <- AGP_unbalanced_rare_high$sample_metadata
count_mat_rare_high <- AGP_unbalanced_rare_high$count_mat
deseq_rare_high <- gmpr_deseq(count_mat_rare_high, sample_metadata_rare_high)
summary_rare_high <- evaluation(taxa_info_unbalanced_rare_high, deseq_rare_high)
saveRDS(summary_rare_high,
        file.path(deseq2_folder, "DA_rare", "high",
                  sprintf("summary_unbalanced_rare_high_%d.rds", ID)))

# abundant DA taxa, low proportion(5%)
AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
                                                 sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
sample_metadata_abundant_low <- AGP_unbalanced_abundant_low$sample_metadata
count_mat_abundant_low <- AGP_unbalanced_abundant_low$count_mat
deseq_abundant_low <- gmpr_deseq(count_mat_abundant_low, sample_metadata_abundant_low)
summary_abundant_low <- evaluation(taxa_info_unbalanced_abundant_low, deseq_abundant_low)
saveRDS(summary_abundant_low,
        file.path(deseq2_folder, "DA_abundant", "low",
                  sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))


# abundant DA taxa, medium proportion(10%)
AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
                                                    sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
sample_metadata_abundant_medium <- AGP_unbalanced_abundant_medium$sample_metadata
count_mat_abundant_medium <- AGP_unbalanced_abundant_medium$count_mat
deseq_abundant_medium <- gmpr_deseq(count_mat_abundant_medium, sample_metadata_abundant_medium)
summary_abundant_medium <- evaluation(taxa_info_unbalanced_abundant_medium, deseq_abundant_medium)
saveRDS(summary_abundant_medium,
        file.path(deseq2_folder, "DA_abundant", "medium",
                  sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))


# abundant DA taxa, high proportion(20%)
AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
                                                  sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
sample_metadata_abundant_high <- AGP_unbalanced_abundant_high$sample_metadata
count_mat_abundant_high <- AGP_unbalanced_abundant_high$count_mat
deseq_abundant_high <- gmpr_deseq(count_mat_abundant_high, sample_metadata_abundant_high)
summary_abundant_high <- evaluation(taxa_info_unbalanced_abundant_high, deseq_abundant_high)
saveRDS(summary_abundant_high,
        file.path(deseq2_folder, "DA_abundant", "high",
                  sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))




