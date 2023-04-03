library(ggplot2)
rm(list=ls())
source('ANCOM+/overdisperse_GLM.R')
source('ANCOM+/Generalized_Least_Squares.R')

otu_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_otu.csv",
                   row.names = 1)

metadata_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_metadata.csv",
                        row.names=1)


glm_result <- pairwise_GLM(count_data=otu_df,
                           metadata=metadata_df,
                           covar="Crohn")

lmer_model <- lmer_function(count_data=otu_df,
                            glm_result = glm_result,
                            reference = "median")

DA_result <- lmer_model[[1]]

write.csv(DA_result, "real_data/Crohn_Vandeputte/Crohn_DiffRatio_result.csv",
          row.names=FALSE)


