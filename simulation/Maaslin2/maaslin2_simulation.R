rm(list=ls())
library(Maaslin2)
library(phyloseq)


data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
maaslin2_folder <- '/home/wangmk/MDAWG/POLDA/simulation/Maaslin2'

source(file.path(maaslin2_folder, "maaslin2_utils.R"))

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ID <- 9

# no DA taxa
AGP_null <- readRDS(file.path(data_folder, "null",
                              sprintf("AGP_simulation_null_%d.rds", ID)))
taxa_info_null <- AGP_null$taxa_info
sample_metadata_null <- AGP_null$sample_metadata
count_mat_null <- AGP_null$count_mat
count_mat_null_df <- as.data.frame(t(count_mat_null))

maaslin2_null <- Maaslin2(input_data=count_mat_null_df,
                           input_metadata = sample_metadata_null,
                           output="tmp",
                           min_prevalence=0.05,
                           normalization="TSS",
                           transform="LOG",
                           analysis_method="LM",
                           fixed_effects="X",
                           standardize = F,
                           cores=2,
                           plot_heatmap = F,
                           plot_scatter = F,
                           correction="BH",
                           max_significance = 0.05)

summary_null <- evaluation(taxa_info_null, maaslin2_null, nullcase=TRUE)
saveRDS(summary_null,
        file.path(maaslin2_folder, "null",
                  sprintf("summary_null_%d.rds", ID)))


# rare DA taxa, low proportion(5%)
AGP_unbalanced_rare_low <- readRDS(file.path(data_folder, "DA_rare", "low",
                                             sprintf("AGP_simulation_unbalanced_rare_low_%d.rds", ID)))
taxa_info_unbalanced_rare_low <- AGP_unbalanced_rare_low$taxa_info
sample_metadata_rare_low <- AGP_unbalanced_rare_low$sample_metadata
count_mat_rare_low <- AGP_unbalanced_rare_low$count_mat
count_mat_rare_low_df <- as.data.frame(t(count_mat_rare_low))

maaslin2_rare_low <- Maaslin2(input_data=count_mat_rare_low_df,
                          input_metadata = sample_metadata_rare_low,
                          output="tmp",
                          min_prevalence=0.05,
                          normalization="TSS",
                          transform="LOG",
                          analysis_method="LM",
                          fixed_effects="X",
                          standardize = F,
                          cores=2,
                          plot_heatmap = F,
                          plot_scatter = F,
                          correction="BH",
                          max_significance = 0.05)

summary_rare_low <- evaluation(taxa_info_unbalanced_rare_low, maaslin2_rare_low)
saveRDS(summary_rare_low,
        file.path(maaslin2_folder, "DA_rare", "low",
                  sprintf("summary_unbalanced_rare_low_%d.rds", ID)))


# rare DA taxa, medium proportion(10%)
AGP_unbalanced_rare_medium <- readRDS(file.path(data_folder, "DA_rare", "medium",
                                             sprintf("AGP_simulation_unbalanced_rare_medium_%d.rds", ID)))
taxa_info_unbalanced_rare_medium <- AGP_unbalanced_rare_medium$taxa_info
sample_metadata_rare_medium <- AGP_unbalanced_rare_medium$sample_metadata
count_mat_rare_medium <- AGP_unbalanced_rare_medium$count_mat
count_mat_rare_medium_df <- as.data.frame(t(count_mat_rare_medium))

maaslin2_rare_medium <- Maaslin2(input_data=count_mat_rare_medium_df,
                              input_metadata = sample_metadata_rare_medium,
                              output="tmp",
                              min_prevalence=0.05,
                              normalization="TSS",
                              transform="LOG",
                              analysis_method="LM",
                              fixed_effects="X",
                              standardize = F,
                              cores=2,
                              plot_heatmap = F,
                              plot_scatter = F,
                              correction="BH",
                              max_significance = 0.05)

summary_rare_medium <- evaluation(taxa_info_unbalanced_rare_medium, maaslin2_rare_medium)
saveRDS(summary_rare_medium,
        file.path(maaslin2_folder, "DA_rare", "medium",
                  sprintf("summary_unbalanced_rare_medium_%d.rds", ID)))


# rare DA taxa, high proportion(20%)
AGP_unbalanced_rare_high <- readRDS(file.path(data_folder, "DA_rare", "high",
                                              sprintf("AGP_simulation_unbalanced_rare_high_%d.rds", ID)))
taxa_info_unbalanced_rare_high <- AGP_unbalanced_rare_high$taxa_info
sample_metadata_rare_high <- AGP_unbalanced_rare_high$sample_metadata
count_mat_rare_high <- AGP_unbalanced_rare_high$count_mat
count_mat_rare_high_df <- as.data.frame(t(count_mat_rare_high))

maaslin2_rare_high <- Maaslin2(input_data=count_mat_rare_high_df,
                               input_metadata = sample_metadata_rare_high,
                               output="tmp",
                               min_prevalence=0.05,
                               normalization="TSS",
                               transform="LOG",
                               analysis_method="LM",
                               fixed_effects="X",
                               standardize = F,
                               cores=2,
                               plot_heatmap = F,
                               plot_scatter = F,
                               correction="BH",
                               max_significance = 0.05)

summary_rare_high <- evaluation(taxa_info_unbalanced_rare_high, maaslin2_rare_high)
saveRDS(summary_rare_high,
        file.path(maaslin2_folder, "DA_rare", "high",
                  sprintf("summary_unbalanced_rare_high_%d.rds", ID)))


# abundant DA taxa, low proportion(5%)
AGP_unbalanced_abundant_low <- readRDS(file.path(data_folder, "DA_abundant", "low",
                                                 sprintf("AGP_simulation_unbalanced_abundant_low_%d.rds", ID)))
taxa_info_unbalanced_abundant_low <- AGP_unbalanced_abundant_low$taxa_info
sample_metadata_abundant_low <- AGP_unbalanced_abundant_low$sample_metadata
count_mat_abundant_low <- AGP_unbalanced_abundant_low$count_mat
count_mat_abundant_low_df <- as.data.frame(t(count_mat_abundant_low))

maaslin2_abundant_low <- Maaslin2(input_data=count_mat_abundant_low_df,
                                  input_metadata = sample_metadata_abundant_low,
                                  output="tmp",
                                  min_prevalence=0.05,
                                  normalization="TSS",
                                  transform="LOG",
                                  analysis_method="LM",
                                  fixed_effects="X",
                                  standardize = F,
                                  cores=2,
                                  plot_heatmap = F,
                                  plot_scatter = F,
                                  correction="BH",
                                  max_significance = 0.05)

summary_abundant_low <- evaluation(taxa_info_unbalanced_abundant_low, maaslin2_abundant_low)
saveRDS(summary_abundant_low,
        file.path(maaslin2_folder, "DA_abundant", "low",
                  sprintf("summary_unbalanced_abundant_low_%d.rds", ID)))


# abundant DA taxa, medium proportion(10%)
AGP_unbalanced_abundant_medium <- readRDS(file.path(data_folder, "DA_abundant", "medium",
                                                    sprintf("AGP_simulation_unbalanced_abundant_medium_%d.rds", ID)))
taxa_info_unbalanced_abundant_medium <- AGP_unbalanced_abundant_medium$taxa_info
sample_metadata_abundant_medium <- AGP_unbalanced_abundant_medium$sample_metadata
count_mat_abundant_medium <- AGP_unbalanced_abundant_medium$count_mat
count_mat_abundant_medium_df <- as.data.frame(t(count_mat_abundant_medium))

maaslin2_abundant_medium <- Maaslin2(input_data=count_mat_abundant_medium_df,
                                     input_metadata = sample_metadata_abundant_medium,
                                     output="tmp",
                                     min_prevalence=0.05,
                                     normalization="TSS",
                                     transform="LOG",
                                     analysis_method="LM",
                                     fixed_effects="X",
                                     standardize = F,
                                     cores=2,
                                     plot_heatmap = F,
                                     plot_scatter = F,
                                     correction="BH",
                                     max_significance = 0.05)

summary_abundant_medium <- evaluation(taxa_info_unbalanced_abundant_medium, maaslin2_abundant_medium)
saveRDS(summary_abundant_medium,
        file.path(maaslin2_folder, "DA_abundant", "medium",
                  sprintf("summary_unbalanced_abundant_medium_%d.rds", ID)))


# abundant DA taxa, high proportion(20%)
AGP_unbalanced_abundant_high <- readRDS(file.path(data_folder, "DA_abundant", "high",
                                                  sprintf("AGP_simulation_unbalanced_abundant_high_%d.rds", ID)))
taxa_info_unbalanced_abundant_high <- AGP_unbalanced_abundant_high$taxa_info
sample_metadata_abundant_high <- AGP_unbalanced_abundant_high$sample_metadata
count_mat_abundant_high <- AGP_unbalanced_abundant_high$count_mat
count_mat_abundant_high_df <- as.data.frame(t(count_mat_abundant_high))

maaslin2_abundant_high <- Maaslin2(input_data=count_mat_abundant_high_df,
                                   input_metadata = sample_metadata_abundant_high,
                                   output="tmp",
                                   min_prevalence=0.05,
                                   normalization="TSS",
                                   transform="LOG",
                                   analysis_method="LM",
                                   fixed_effects="X",
                                   standardize = F,
                                   cores=2,
                                   plot_heatmap = F,
                                   plot_scatter = F,
                                   correction="BH",
                                   max_significance = 0.05)

summary_abundant_high <- evaluation(taxa_info_unbalanced_abundant_high, maaslin2_abundant_high)
saveRDS(summary_abundant_high,
        file.path(maaslin2_folder, "DA_abundant", "high",
                  sprintf("summary_unbalanced_abundant_high_%d.rds", ID)))


