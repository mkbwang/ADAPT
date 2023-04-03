

rm(list=ls())
source('ANCOM+/overdisperse_GLM.R')
source('ANCOM+/Generalized_Least_Squares.R')

otu_infant <- read.csv("real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_infant_OTU.csv",
                   row.names = 1)

metadata_infant <- read.csv("real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_infant_metadata.csv",
                        row.names=1)
metadata_infant$Malawi <- metadata_infant$country == "MA"
glm_result_infant <- pairwise_GLM(count_data=otu_infant,
                           meta=metadata_infant,
                           covar="Malawi")

lmer_model <- lmer_function(count_data=otu_infant,
                         glm_result = glm_result_infant,
                         reference = "median")

DA_result_infant <- lmer_model[[1]]
write.csv(DA_result_infant, "real_data/Yatsunenko_data/Malawi_Venezuela/MA_VEN_infant_DiffRatio_result.csv",
          row.names = FALSE)

