rm(list=ls())

input_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
output_folder <- '/home/wangmk/MDAWG/POLDA/simulation/POLDA'

source("/home/wangmk/MDAWG/POLDA/POLDA/POLDA.R")

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(input_folder, input_fname))

count_mat <- simulated_data$otu_tab_sim
covariate <- simulated_data$metadata
truth <- simulated_data$taxa

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

