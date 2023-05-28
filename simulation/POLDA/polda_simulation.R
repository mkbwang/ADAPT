rm(list=ls())
library(phyloseq)
library(ClassComparison)
library(DescTools)
library(qvalue)
library(vegan)
library(ALDEx2)
library(dirmult)
input_folder <- 'simulation/data/'
output_folder <- 'simulation/POLDA'

source("POLDA/POLDA.R")

template <- readRDS(file.path(input_folder, 'AGP_template.rds'))

original_count_mat <- otu_table(template)@.Data
prevalence <- rowMeans(original_count_mat != 0)
metadata <- sample_data(template) |> data.frame()

# shuffle the counts within a taxon
permuted_count <- original_count_mat
Ntaxa <- nrow(original_count_mat)
Nindv <- ncol(original_count_mat)
for (j in 1:Ntaxa){
  permuted_count[j, ] <- original_count_mat[j, sample(1:Nindv, Nindv)]
}
permuted_depths <- colSums(permuted_count)
hist(permuted_depths, nclass=20)


# result_original <- polda(otu_table = original_count_mat,
#                          metadata = metadata,
#                          covar="X")

relabd_result_permute <- reference_GLM(permuted_count,
                                    metadata=metadata,
                                    reftaxa = rownames(permuted_count),
                                    covar="X")
hist(relabd_result_permute$pval, nclass=30)
pvals <- relabd_result_permute$pval

bumfit <- Bum(pvals)
lambda_hat <- bumfit@lhat
a_hat <- bumfit@ahat

# likelihood of BUM
llk_func <- function(lambda_hat, a_hat, p_val){
  result <- lambda_hat + (1-lambda_hat) * a_hat * (p_val)^(a_hat - 1)
  return(result)
}
# Cutoff calculation for FDR
cutoff_func <- function(lambda_hat, a_hat, alpha){
  pi_hat <- lambda_hat + (1 - lambda_hat) * a_hat
  result <- ((pi_hat - alpha*lambda_hat)/(alpha*(1-lambda_hat)))^(1/(a_hat - 1))
  return(result)
}

llk_vals <- rep(0, 1000)
for (j in 1:1000){
  llk_vals[j] <- llk_func(lambda_hat=bumfit@lhat,
                          a_hat=bumfit@ahat,
                          p_val=pvals[j])
}
package_llk <- likelihoodBum(bumfit)
stopifnot(all(llk_vals == package_llk))

manual_cutoff <- cutoff_func(lambda_hat, a_hat, 0.05)
package_cutoff <- cutoffSignificant(bumfit, alpha=0.05, by="FDR")
stopifnot(manual_cutoff == package_cutoff)


result_permuted <- polda(otu_table = permuted_count,
                         metadata=metadata,
                         covar="X")



# ID <- 50
# input_fname <- sprintf('simulated_%d.rds', ID)
# simulated_data <- readRDS(file.path(input_folder, input_fname))
# count_mat <- simulated_data$count_mat
# metadata <- simulated_data$metadata
# taxa_info <- simulated_data$taxa_info
# result <- polda(count_mat, metadata, 'group', covartype="categorical")
# 
# phy_obj <- phyloseq(otu_table(count_mat, taxa_are_rows=TRUE),
#                     sample_data(metadata))
# 
# ancombc_result <- ancombc(phyloseq=phy_obj, formula="group",
#                           p_adj_method = "BH", zero_cut=0.95, group="group")
# aldex_result <- aldex(count_mat, metadata$group, test="t")
# 
# taxa_info_df <- as.data.frame(taxa_info)
# taxa_info_df$polda <- FALSE
# taxa_info_df$ancombc <- FALSE
# taxa_info_df[result$DA_taxa, 'polda'] <- TRUE
# taxa_info_df$ancombc <- ancombc_result$res$diff_abn$group
# 
# table(taxa_info_df$fold_change_e == 1, taxa_info_df$polda)
# table(taxa_info_df$fold_change_e == 1, taxa_info_df$ancombc)

for (ID in 1:100){
  print(ID)
  input_fname <- sprintf('simulated_data_%d.rds', ID)
  simulated_data <- readRDS(file.path(input_folder, input_fname))
  
  count_mat <- simulated_data$otu_tab_sim
  covariate <- simulated_data$metadata
  truth <- simulated_data$taxa
  
  
  result <- polda(count_mat, covariate, 'X', covartype="categorical")
  
  #truth_subset <- truth %>% filter(!Taxa %in% reftaxa)
  truth$Taxon <- rownames(truth)
  pval_df <- result$P_Value
  performance <- truth %>% left_join(pval_df, by="Taxon")
  performance$RefTaxa <- FALSE
  performance$RefTaxa[performance$Taxon %in% result$Reference_Taxa] <- TRUE
  
  write.csv(performance, 
            file.path(output_folder, sprintf('polda_simulation_result_%d.csv', ID)),
            row.names=FALSE)
}


