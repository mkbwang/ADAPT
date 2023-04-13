library(ggplot2)
library(dplyr)
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

# glm_result <- read.csv("real_data/Crohn_Vandeputte/Crohn_GLM_perm0.csv")

full_model <- generalized_ls(otu_df, glm_result)
full_model_params <- full_model$fitted_parameters
write.csv(full_model_params,
          'real_data/Crohn_Vandeputte/Crohn_DiffRatio_result.csv', row.names=FALSE)
tau_abs <- abs(full_model$fitted_parameters$tau_hat)
teststat_abs <- abs(full_model$fitted_parameters$teststat)
backward_selection_order <- order(teststat_abs)
backward_selection_order2 <- order(tau_abs)
taxa_names <- row.names(otu_df)


AIC_scores <- rep(0, 30)
AIC_scores[1] <- full_model$AIC_score
small_models <- list()
for (j in 1:29){
  small_models[[j]] <- generalized_ls(otu_df, glm_result,
                                      shrinkid = backward_selection_order[seq(1, j)])
  AIC_scores[j+1] <- small_models[[j]]$AIC_score
}


ref_taxa <- taxa_names[backward_selection_order[seq(1, which.min(AIC_scores) - 1)]]

# select subsets of
glm_filter <- xor(glm_result$T1 %in% ref_taxa, glm_result$T2 %in% ref_taxa)
glmresult_subset <- glm_result[glm_filter, ]

glmresult_subset$pval_adjust <- p.adjust(glmresult_subset$pval,
                                         method="BH")

glmresult_subset_significant <- glmresult_subset[glmresult_subset$pval_adjust < 0.05, ] |>
  na.omit()

difftaxa <- setdiff(c(glmresult_subset_significant$T1, glmresult_subset_significant$T2),
                    ref_taxa)

result <- list(reftaxa = ref_taxa,
               difftaxa = difftaxa)

saveRDS(result, "real_data/Crohn_Vandeputte/Crohn_DiffRatio_result.rds")

