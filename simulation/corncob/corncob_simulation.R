rm(list=ls())
library(corncob)
library(phyloseq)


data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
corncob_folder <- '/home/wangmk/MDAWG/POLDA/simulation/corncob'
source(file.path(corncob_folder, 'corncob_utils.R'))

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 10

# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata_null <- AGP_null$sample_metadata
phyobj_null <- phyloseq(otu_table(AGP_null$count_mat, taxa_are_rows=TRUE),
                        sample_data(sample_metadata_null))
corncob_null <- differentialTest(formula=~X,
                                     phi.formula=~1,
                                     formula_null=~1,
                                     phi.formula_null=~1,
                                     data=phyobj_null,
                                     fdr_cutoff=0.05,
                                     test="Wald")
summary_null <- evaluation(taxa_info_null, corncob_null, nullcase=TRUE)
saveRDS(summary_null,
        file.path(corncob_folder, "null",
                  sprintf("summary_null_%d.rds", ID)))


# rare DA taxa, low proportion(5%)
AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
                                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
sample_metadata_rare_low <- AGP_unbalanced_rare_low$sample_metadata
phyobj_unbalanced_rare_low <- phyloseq(otu_table(AGP_unbalanced_rare_low$count_mat, taxa_are_rows = TRUE),
                                       sample_data(sample_metadata_rare_low))
corncob_rare_low <- differentialTest(formula=~X,
                                     phi.formula=~1,
                                     formula_null=~1,
                                     phi.formula_null = ~1,
                                     data=phyobj_unbalanced_rare_low,
                                     fdr_cutoff=0.05,
                                     test="Wald")
summary_rare_low <- evaluation(taxa_info_unbalanced_rare_low, 
                               corncob_rare_low)
saveRDS(summary_rare_low,
        file.path(corncob_folder, "DA_rare", "low",
                  sprintf("summary_unbalanced_rare_low_%d.rds", ID)))



# rare DA taxa, medium proportion(10%)
AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                             sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
sample_metadata_rare_medium <- AGP_unbalanced_rare_medium$sample_metadata
phyobj_unbalanced_rare_medium <- phyloseq(otu_table(AGP_unbalanced_rare_medium$count_mat, taxa_are_rows = TRUE),
                                       sample_data(sample_metadata_rare_medium))
corncob_rare_medium <- differentialTest(formula=~X,
                                     phi.formula=~1,
                                     formula_null=~1,
                                     phi.formula_null = ~1,
                                     data=phyobj_unbalanced_rare_medium,
                                     fdr_cutoff=0.05,
                                     test="Wald")
summary_rare_medium <- evaluation(taxa_info_unbalanced_rare_medium, 
                               corncob_rare_medium)
saveRDS(summary_rare_medium,
        file.path(corncob_folder, "DA_rare", "medium",
                  sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))


# rare DA taxa, high proportion(20%)
AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                                sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
sample_metadata_rare_high <- AGP_unbalanced_rare_high$sample_metadata
phyobj_unbalanced_rare_high <- phyloseq(otu_table(AGP_unbalanced_rare_high$count_mat, taxa_are_rows = TRUE),
                                          sample_data(sample_metadata_rare_high))
corncob_rare_high <- differentialTest(formula=~X,
                                        phi.formula=~1,
                                        formula_null=~1,
                                        phi.formula_null = ~1,
                                        data=phyobj_unbalanced_rare_high,
                                        fdr_cutoff=0.05,
                                        test="Wald")
summary_rare_high <- evaluation(taxa_info_unbalanced_rare_high, 
                                  corncob_rare_high)
saveRDS(summary_rare_high,
        file.path(corncob_folder, "DA_rare", "high",
                  sprintf("summary_unbalanced_rare_high_%d.rds", ID)))


# abundant DA taxa, low proportion(5%)
AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
                                                 sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
sample_metadata_abundant_low <- AGP_unbalanced_abundant_low$sample_metadata
phyobj_unbalanced_abundant_low <- phyloseq(otu_table(AGP_unbalanced_abundant_low$count_mat, taxa_are_rows = TRUE),
                                           sample_data(sample_metadata_abundant_low))
corncob_abundant_low <- differentialTest(formula=~X,
                                         phi.formula=~1,
                                         formula_null=~1,
                                         phi.formula_null = ~1,
                                         data=phyobj_unbalanced_abundant_low,
                                         fdr_cutoff=0.05,
                                         test="Wald")
summary_abundant_low <- evaluation(taxa_info_unbalanced_abundant_low, 
                                   corncob_abundant_low)
saveRDS(summary_abundant_low,
        file.path(corncob_folder, "DA_abundant", "low",
                  sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))


# abundant DA taxa, medium proportion(10%)
AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
                                                    sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
sample_metadata_abundant_medium <- AGP_unbalanced_abundant_medium$sample_metadata
phyobj_unbalanced_abundant_medium <- phyloseq(otu_table(AGP_unbalanced_abundant_medium$count_mat, taxa_are_rows = TRUE),
                                              sample_data(sample_metadata_abundant_medium))
corncob_abundant_medium <- differentialTest(formula=~X,
                                            phi.formula=~1,
                                            formula_null=~1,
                                            phi.formula_null = ~1,
                                            data=phyobj_unbalanced_abundant_medium,
                                            fdr_cutoff=0.05,
                                            test="Wald")
summary_abundant_medium <- evaluation(taxa_info_unbalanced_abundant_medium, 
                                      corncob_abundant_medium)
saveRDS(summary_abundant_medium,
        file.path(corncob_folder, "DA_abundant", "medium",
                  sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))


# abundant DA taxa, high proportion(20%)
AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
                                                  sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
sample_metadata_abundant_high <- AGP_unbalanced_abundant_high$sample_metadata
phyobj_unbalanced_abundant_high <- phyloseq(otu_table(AGP_unbalanced_abundant_high$count_mat, taxa_are_rows = TRUE),
                                            sample_data(sample_metadata_abundant_high))
corncob_abundant_high <- differentialTest(formula=~X,
                                          phi.formula=~1,
                                          formula_null=~1,
                                          phi.formula_null = ~1,
                                          data=phyobj_unbalanced_abundant_high,
                                          fdr_cutoff=0.05,
                                          test="Wald")
summary_abundant_high <- evaluation(taxa_info_unbalanced_abundant_high, 
                                    corncob_abundant_high)
saveRDS(summary_abundant_high,
        file.path(corncob_folder, "DA_abundant", "high",
                  sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))

