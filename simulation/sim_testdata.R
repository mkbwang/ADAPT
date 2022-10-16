rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'


num_taxa = 1000
sample_num_grp1 <- 50
sample_num_grp2 <- 50
low_abundance <- 50
mid_abundance <- 200
high_abundance <- 10000
proportion_difference <- 0.25
seed <- 1
struc_zeroprop <- 0
outlier_zeroprop <- 0
balanced_microbiome_load <- FALSE
balanced_libsize <- TRUE
sample_fraction <- 'small'

library(phyloseq)
source(file.path(folder, 'sim_data_poi_gam_two_grp.R'))


for (seed in c(47)){
  print(seed)
  simulated_data <- abn.tab.gen1(num_taxa, sample_num_grp1, sample_num_grp2,
                                 low_abundance, mid_abundance, high_abundance,
                                 proportion_difference, seed, seed, struc_zeroprop,
                                 outlier_zeroprop, balanced_microbiome_load, balanced_libsize,
                                 sample_fraction)
  # abn_summary <- simulated_data$mean.eco.abn
  # obs_abundance <- simulated_data$obs.abn
  # abn_summary$effect.size <- abn_summary$temp.grp1 / abn_summary$temp.grp2
  # simulated_data$mean.eco.abn <- abn_summary

  filename <- sprintf('simulated_data_%d.rds', seed)
  saveRDS(simulated_data, file=file.path(folder, 'data', filename))
}


# take a small subset of taxa who are abundant
# selected_taxa <- which(abn_summary$temp.grp1 > 200 & abn_summary$temp.grp2 > 200)
#
# subset_abn_summary <- abn_summary[selected_taxa, ]
# subset_sampled_count <- simulated_data$obs.abn[selected_taxa, ]
# indv_group <- simulated_data$grp
#
# subset_data <- list(abn_info = subset_abn_summary,
#                     sample_counts = subset_sampled_count,
#                     indv_metadata = indv_group)
# saveRDS(subset_data, file=file.path(folder, 'subset_simulated_data_1.rds'))


