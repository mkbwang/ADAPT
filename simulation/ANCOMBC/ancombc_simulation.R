rm(list=ls())
library(ANCOMBC)
library(phyloseq)

data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
ancombc_folder <- '/home/wangmk/MDAWG/POLDA/simulation/ANCOMBC'

source(file.path(ancombc_folder, "ancombc_utils.R"))
ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 10

# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata_null <- AGP_null$sample_metadata
phyobj_null <- phyloseq(otu_table(AGP_null$count_mat, taxa_are_rows=TRUE),
                        sample_data(sample_metadata_null))
ancombc_null <- ancombc(phyobj_null, formula = 'X', p_adj_method='BH',
                          group='X', struc_zero=TRUE, prv_cut=0.05)
performance_null <- evaluation(taxa_truth=taxa_info_null,
                               ancombc_result=ancombc_null,
                               nullcase=TRUE)
output_null <- list(performance=performance_null,
                    ancombc_result=ancombc_null)
saveRDS(output_null,
        file.path(ancombc_folder, "null",
                  sprintf("summary_null_%d.rds", ID)))


# rare DA taxa, low proportion(5%)
AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
                                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
sample_metadata_rare_low <- AGP_unbalanced_rare_low$sample_metadata
phyobj_unbalanced_rare_low <- phyloseq(otu_table(AGP_unbalanced_rare_low$count_mat, taxa_are_rows = TRUE),
                                       sample_data(sample_metadata_rare_low))
ancombc_unbalanced_rare_low <- ancombc(phyobj_unbalanced_rare_low, formula = 'X', p_adj_method='BH',
                                       group='X', struc_zero=TRUE, prv_cut=0.05)
performance_unbalanced_rare_low <- evaluation(taxa_info_unbalanced_rare_low,
                                              ancombc_unbalanced_rare_low)
output_unbalanced_rare_low <- list(performance=performance_unbalanced_rare_low,
                                   ancombc_result=ancombc_unbalanced_rare_low)
saveRDS(output_unbalanced_rare_low,
        file.path(ancombc_folder, "DA_rare", "low",
                  sprintf("summary_unbalanced_rare_low_%d.rds", ID)))


# rare DA taxa, medium proportion(10%)
AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                             sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
sample_metadata_rare_medium <- AGP_unbalanced_rare_medium$sample_metadata
phyobj_unbalanced_rare_medium <- phyloseq(otu_table(AGP_unbalanced_rare_medium$count_mat, taxa_are_rows = TRUE),
                                       sample_data(sample_metadata_rare_medium))
ancombc_unbalanced_rare_medium <- ancombc(phyobj_unbalanced_rare_medium, formula = 'X', p_adj_method='BH',
                                       group='X', struc_zero=TRUE, prv_cut=0.05)
performance_unbalanced_rare_medium <- evaluation(taxa_info_unbalanced_rare_medium,
                                              ancombc_unbalanced_rare_medium)
output_unbalanced_rare_medium <- list(performance=performance_unbalanced_rare_medium,
                                   ancombc_result=ancombc_unbalanced_rare_medium)
saveRDS(output_unbalanced_rare_medium,
        file.path(ancombc_folder, "DA_rare", "medium",
                  sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))


# rare DA taxa, high proportion (20%)
AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                                sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
sample_metadata_rare_high <- AGP_unbalanced_rare_high$sample_metadata
phyobj_unbalanced_rare_high <- phyloseq(otu_table(AGP_unbalanced_rare_high$count_mat, taxa_are_rows = TRUE),
                                          sample_data(sample_metadata_rare_high))
ancombc_unbalanced_rare_high <- ancombc(phyobj_unbalanced_rare_high, formula = 'X', p_adj_method='BH',
                                          group='X', struc_zero=TRUE, prv_cut=0.05)
performance_unbalanced_rare_high <- evaluation(taxa_info_unbalanced_rare_high,
                                                 ancombc_unbalanced_rare_high)
output_unbalanced_rare_high <- list(performance=performance_unbalanced_rare_high,
                                      ancombc_result=ancombc_unbalanced_rare_high)
saveRDS(output_unbalanced_rare_high,
        file.path(ancombc_folder, "DA_rare", "high",
                  sprintf("summary_unbalanced_rare_high_%d.rds", ID)))



# abundant DA taxa, low proportion (5%)
AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
                                              sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
sample_metadata_abundant_low <- AGP_unbalanced_abundant_low$sample_metadata
phyobj_unbalanced_abundant_low <- phyloseq(otu_table(AGP_unbalanced_abundant_low$count_mat, taxa_are_rows = TRUE),
                                        sample_data(sample_metadata_abundant_low))
ancombc_unbalanced_abundant_low <- ancombc(phyobj_unbalanced_abundant_low, formula = 'X', p_adj_method='BH',
                                        group='X', struc_zero=TRUE, prv_cut=0.05)
performance_unbalanced_abundant_low <- evaluation(taxa_info_unbalanced_abundant_low,
                                               ancombc_unbalanced_abundant_low)
output_unbalanced_abundant_low <- list(performance=performance_unbalanced_abundant_low,
                                    ancombc_result=ancombc_unbalanced_abundant_low)
saveRDS(output_unbalanced_abundant_low,
        file.path(ancombc_folder, "DA_abundant", "low",
                  sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))


# abundant DA taxa, medium proportion (10%)
AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
                                                    sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
sample_metadata_abundant_medium <- AGP_unbalanced_abundant_medium$sample_metadata
phyobj_unbalanced_abundant_medium <- phyloseq(otu_table(AGP_unbalanced_abundant_medium$count_mat, taxa_are_rows = TRUE),
                                              sample_data(sample_metadata_abundant_medium))
ancombc_unbalanced_abundant_medium <- ancombc(phyobj_unbalanced_abundant_medium, formula = 'X', p_adj_method='BH',
                                              group='X', struc_zero=TRUE, prv_cut=0.05)
performance_unbalanced_abundant_medium <- evaluation(taxa_info_unbalanced_abundant_medium,
                                                     ancombc_unbalanced_abundant_medium)
output_unbalanced_abundant_medium <- list(performance=performance_unbalanced_abundant_medium,
                                          ancombc_result=ancombc_unbalanced_abundant_medium)
saveRDS(output_unbalanced_abundant_medium,
        file.path(ancombc_folder, "DA_abundant", "medium",
                  sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))


# abundant DA taxa, high proportion (20%)

AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
                                                  sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
sample_metadata_abundant_high <- AGP_unbalanced_abundant_high$sample_metadata
phyobj_unbalanced_abundant_high <- phyloseq(otu_table(AGP_unbalanced_abundant_high$count_mat, taxa_are_rows = TRUE),
                                            sample_data(sample_metadata_abundant_high))
ancombc_unbalanced_abundant_high <- ancombc(phyobj_unbalanced_abundant_high, formula = 'X', p_adj_method='BH',
                                            group='X', struc_zero=TRUE, prv_cut=0.05)
performance_unbalanced_abundant_high <- evaluation(taxa_info_unbalanced_abundant_high,
                                                   ancombc_unbalanced_abundant_high)
output_unbalanced_abundant_high <- list(performance=performance_unbalanced_abundant_high,
                                        ancombc_result=ancombc_unbalanced_abundant_high)
saveRDS(output_unbalanced_abundant_high, 
        file.path(ancombc_folder, "DA_abundant", "high", 
                  sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))

