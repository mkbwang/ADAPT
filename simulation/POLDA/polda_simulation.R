rm(list=ls())
library(phyloseq)
library(BioNet)
library(DescTools)
library(qvalue)
library(vegan)
library(Rfast)
library(ALDEx2)
input_folder <- 'simulation/data/'
output_folder <- 'simulation/POLDA'

source("POLDA/POLDA.R")


template <- readRDS(file.path(input_folder, 'Yatsunenko_template.rds'))

original_count_mat <- otu_table(template)@.Data
prevalence <- rowMeans(original_count_mat != 0)
metadata <- sample_data(template) |> as.data.frame()
metadata$group <- c(rep(0, 64), rep(1, 65))

# shuffle the counts within a taxon
permuted_count <- original_count_mat
Ntaxa <- nrow(original_count_mat)
Nindv <- ncol(original_count_mat)
for (j in 1:Ntaxa){
  permuted_count[j, ] <- original_count_mat[j, sample(1:Nindv, Nindv)]
}

relabd_result_permute <- relabd_GLM(permuted_count,
                                    metadata=metadata,
                                    covar="group")
hist(relabd_result_permute$pval, nclass=50)
AndersonDarlingTest(relabd_result_permute$pval, null="punif")
ks.test(relabd_result_permute$pval, y="punif")

BH_pval <- p.adjust(relabd_result_permute$pval, method="BH")
qvalues <- qvalue(relabd_result_permute$pval)

median_effect <- median(relabd_result$effect)
effect_rank <- rank(abs(relabd_result$effect - median(relabd_result$effect)))

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


