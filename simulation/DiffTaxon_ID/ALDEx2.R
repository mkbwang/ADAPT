rm(list=ls())
library(dplyr)
library(ALDEx2)


folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'

ID <- 1
input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(folder, 'data', 'semiparametric', input_fname))

count_mat <- simulated_data$otu.tab.sim
covariate <- as.vector(simulated_data$covariate)

aldex_result <- aldex(count_mat, covariate, test='t')

result_df <- data.frame(DA=aldex_result$wi.eBH < 0.05,
                        truth = simulated_data$diff.otu.ind)

row.names(result_df) <- row.names(count_mat)

