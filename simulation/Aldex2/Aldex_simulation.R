rm(list=ls())
library(dplyr)
library(ALDEx2)

input_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
output_folder <- '/home/wangmk/MDAWG/POLDA/simulation/Aldex2'

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(input_folder, input_fname))

count_mat <- simulated_data$otu_tab_sim
covariate <- as.vector(simulated_data$metadata$X)
truth <- simulated_data$taxa

aldex_result <- aldex(count_mat, covariate, test='t')

performance <- data.frame(Taxon=rownames(truth),
                        DA_Pval=aldex_result$wi.eBH,
                        logfold=truth$logfold)

write.csv(performance, 
          file.path(output_folder, sprintf('aldex_simulation_result_%d.csv', ID)),
          row.names=FALSE)

