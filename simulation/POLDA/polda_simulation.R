rm(list=ls())

input_folder <- 'simulation/data/semiparametric'
output_folder <- '/home/wangmk/MDAWG/POLDA/simulation/POLDA'

source("POLDA/POLDA.R")

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ID <- 27
input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(input_folder, input_fname))

count_mat <- simulated_data$otu_tab_sim
covariate <- simulated_data$metadata
truth <- simulated_data$taxa

alltaxa_result <- relabd_GLM(count_data=count_mat,
                             metadata=covariate,
                             covar="X")
alltaxa_result$pval_adjust <- p.adjust(alltaxa_result$pval, method="BH")

med_relabd_odds <- median(alltaxa_result$effect)

# first check the "depleted" taxa group
Taxa_D <- alltaxa_result$Taxa[alltaxa_result$effect < med_relabd_odds]
while(1) {
  count_mat_D <- count_mat[Taxa_D, ]
  deplete_taxa_result <- relabd_GLM(count_data=count_mat_D,
                                    metadata=covariate,
                                    covar="X")
  if (sum(deplete_taxa_result$teststat < qnorm(0.05)) == 0) break
  # if there are traces of DA taxa, cut away those whose estimated beta are smaller than zero
  removed_taxa <- deplete_taxa_result$Taxa[deplete_taxa_result$effect < 0]
  print(removed_taxa)
  Taxa_D <- setdiff(Taxa_D, removed_taxa)
}


Taxa_E <- alltaxa_result$Taxa[alltaxa_result$effect > med_relabd_odds]
while(1) {
  count_mat_E <- count_mat[Taxa_E, ]
  enriched_taxa_result <- relabd_GLM(count_data=count_mat_E,
                                    metadata=covariate,
                                    covar="X")
  if (sum(enriched_taxa_result$teststat > -qnorm(0.05)) == 0) break
  removed_taxa <- enriched_taxa_result$Taxa[enriched_taxa_result$effect > 0]
  print(removed_taxa)
  Taxa_E <- setdiff(Taxa_E, removed_taxa)
}

reference_taxa <- c(Taxa_D, Taxa_E)

refglm_result <- reference_GLM(count_data = count_mat, 
                               metadata = covariate,
                               covar="X",
                               reftaxa=reference_taxa)

result <- polda(count_mat, covariate, 'X', covartype="categorical",
                startdrop="median")

#truth_subset <- truth %>% filter(!Taxa %in% reftaxa)
truth$Taxon <- rownames(truth)
pval_df <- result$P_Value
performance <- truth %>% left_join(pval_df, by="Taxon")
performance$RefTaxa <- FALSE
performance$RefTaxa[performance$Taxon %in% result$Reference_Taxa] <- TRUE

write.csv(performance, 
          file.path(output_folder, sprintf('polda_simulation_result_%d.csv', ID)),
          row.names=FALSE)

