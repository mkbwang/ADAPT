rm(list=ls())
library(eBay)

data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
ebay_folder <- '/home/wangmk/MDAWG/POLDA/simulation/eBay'

source(file.path(ebay_folder, "ebay_utils.R"))

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 9


# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata_null <- AGP_null$sample_metadata
count_mat_null <- AGP_null$count_mat

eBay_null <- eBay(otu.data=t(count_mat_null), group=sample_metadata_null$X,
                  test.method="wilcoxon", cutf=0.05, adj.m="none")
summary_null <- evaluation(taxa_info_null, eBay_null, nullcase = TRUE)
saveRDS(summary_null,
        file.path(ebay_folder, "null",
                  sprintf("summary_null_%d.rds", ID)))



# rare DA taxa, low proportion(5%)
AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
                                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
sample_metadata_rare_low <- AGP_unbalanced_rare_low$sample_metadata
count_mat_rare_low <- AGP_unbalanced_rare_low$count_mat
eBay_rare_low <- eBay(otu.data=t(count_mat_rare_low), group=sample_metadata_rare_low$X,
                      test.method="wilcoxon", cutf=0.05, adj.m="none")
summary_rare_low <- evaluation(taxa_info_unbalanced_rare_low, eBay_rare_low)
saveRDS(summary_rare_low,
        file.path(ebay_folder, "DA_rare", "low",
                  sprintf("summary_unbalanced_rare_low_%d.rds", ID)))


# rare DA taxa, medium proportion(10%)
AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                                sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
sample_metadata_rare_medium <- AGP_unbalanced_rare_medium$sample_metadata
count_mat_rare_medium <- AGP_unbalanced_rare_medium$count_mat
eBay_rare_medium <- eBay(otu.data=t(count_mat_rare_medium), group=sample_metadata_rare_medium$X,
                         test.method="wilcoxon", cutf=0.05, adj.m="none")
summary_rare_medium <- evaluation(taxa_info_unbalanced_rare_medium, eBay_rare_medium)
saveRDS(summary_rare_medium,
        file.path(ebay_folder, "DA_rare", "medium",
                  sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))


# rare DA taxa, high proportion(20%)
AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                              sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
sample_metadata_rare_high <- AGP_unbalanced_rare_high$sample_metadata
count_mat_rare_high <- AGP_unbalanced_rare_high$count_mat
eBay_rare_high <- eBay(otu.data=t(count_mat_rare_high), group=sample_metadata_rare_high$X,
                       test.method="wilcoxon", cutf=0.05, adj.m="none")
summary_rare_high <- evaluation(taxa_info_unbalanced_rare_high, eBay_rare_high)
saveRDS(summary_rare_high,
        file.path(ebay_folder, "DA_rare", "high",
                  sprintf("summary_unbalanced_rare_high_%d.rds", ID)))


# abundant DA taxa, low proportion(5%)
AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
                                                 sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
sample_metadata_abundant_low <- AGP_unbalanced_abundant_low$sample_metadata
count_mat_abundant_low <- AGP_unbalanced_abundant_low$count_mat
eBay_abundant_low <- eBay(otu.data=t(count_mat_abundant_low), group=sample_metadata_abundant_low$X,
                          test.method="wilcoxon", cutf=0.05, adj.m="none")
summary_abundant_low <- evaluation(taxa_info_unbalanced_abundant_low, eBay_abundant_low)
saveRDS(summary_abundant_low,
        file.path(ebay_folder, "DA_abundant", "low",
                  sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))


# abundant DA taxa, medium proportion(10%)
AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
                                                    sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
sample_metadata_abundant_medium <- AGP_unbalanced_abundant_medium$sample_metadata
count_mat_abundant_medium <- AGP_unbalanced_abundant_medium$count_mat
eBay_abundant_medium <- eBay(otu.data=t(count_mat_abundant_medium), group=sample_metadata_abundant_medium$X,
                             test.method="wilcoxon", cutf=0.05, adj.m="none")
summary_abundant_medium <- evaluation(taxa_info_unbalanced_abundant_medium, eBay_abundant_medium)
saveRDS(summary_abundant_medium,
        file.path(ebay_folder, "DA_abundant", "medium",
                  sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))


# abundant DA taxa, high proportion(20%)
AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
                                                  sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
sample_metadata_abundant_high <- AGP_unbalanced_abundant_high$sample_metadata
count_mat_abundant_high <- AGP_unbalanced_abundant_high$count_mat
eBay_abundant_high <- eBay(otu.data=t(count_mat_abundant_high), group=sample_metadata_abundant_high$X,
                           test.method="wilcoxon", cutf=0.05, adj.m="none")
summary_abundant_high <- evaluation(taxa_info_unbalanced_abundant_high, eBay_abundant_high)
saveRDS(summary_abundant_high,
        file.path(ebay_folder, "DA_abundant", "high",
                  sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))



