library(ggplot2)
library(dplyr)
library(reshape2)
rm(list=ls())
source('POLDA/POLDA.R')

otu_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_otu.csv",
                   row.names = 1)

metadata_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_metadata.csv",
                        row.names=1)

result <- polda(otu_table=otu_df, metadata=metadata_df, covar="Crohn")

saveRDS(result, "real_data/Crohn_Vandeputte/Crohn_DiffRatio_result.rds")

