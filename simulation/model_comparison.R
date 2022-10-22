rm(list=ls())

folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
library(dplyr)


fnum <- 1
data_name <- sprintf('simulated_data_%d.rds', fnum)
data <- readRDS(file.path(folder, 'data', data_name))
abn_info <- data$mean.eco.abn

taxa_names <- row.names(abn_info)
taxa_pairs <- combn(taxa_names, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("T1", "T2")

# collect the truth
truth_vec <- c()
for (j in 1:999){
  ref_taxon <- taxa_names[j]
  ref_effsize <- abn_info[j, 'effect.size']
  others_effsize <- abn_info[(j+1):1000, 'effect.size']
  truth_vec <- c(truth_vec, others_effsize != ref_effsize)
}


# load model results
ANCOM_file <- sprintf('ancom_result_%d.csv', fnum)
ANCOM_result <- read.csv(file.path(folder, 'ancom_result', ANCOM_file))

ZOIB_file <- sprintf('zoib_result_%d.csv', fnum)
ZOIB_result <- read.csv(file.path(folder, 'ZOIB_result', ZOIB_file))

GLMM_file <- sprintf('glmm_result_%d.csv', fnum)
GLMM_result <- read.csv(file.path(folder, 'glmm_result', GLMM_file))

combination <- cbind(taxa_pairs, truth_vec, ANCOM_result$pval, ZOIB_result$zinfpval,
                     ZOIB_result$oinfpval, ZOIB_result$betapval, GLMM_result$pval)
colnames(combination) <- c("T1", "T2", "diffratio", "ANCOM_pval",
                            'ZOIB_zinf_pval', 'ZOIB_oinf_pval',
                           'ZOIB_beta_pval', 'GLMM_pval')


combination_replicate <- combination

combination_replicate$ANCOM_pval <- p.adjust(combination_replicate$ANCOM_pval,
                                             method='BH')
combination_replicate$GLMM_pval <- p.adjust(combination_replicate$GLMM_pval,
                                            method='BH')
combination_replicate$ZOIB_zinf_pval <- p.adjust(combination_replicate$ZOIB_zinf_pval,
                                            method='BH')
combination_replicate$ZOIB_oinf_pval <- p.adjust(combination_replicate$ZOIB_oinf_pval,
                                                 method='BH')
combination_replicate$ZOIB_beta_pval <- p.adjust(combination_replicate$ZOIB_beta_pval,
                                                 method='BH')



## ANCOM performance
ancom_decision <- combination_replicate$ANCOM_pval <0.05

map_significance <- function(pairresult, num_taxa=1000){
  decision_mat <- matrix(FALSE, nrow=num_taxa, ncol=num_taxa)
  cursor <- 1
  for (j in 1:(num_taxa-1)){
    decision_mat[j, (j+1):num_taxa] <- pairresult[cursor: (cursor+num_taxa-j-1)]
    cursor <- cursor + num_taxa -j
  }
  decision_mat <- decision_mat + t(decision_mat)
  return(decision_mat)
}

# GLMM_decision <- combination_replicate$GLMM_pval <0.05
difftaxa_truth <- abn_info$effect.size != 1

library(pROC)
ancom_decision_mat <- map_significance(ancom_decision)
ancom_taxon_prob <- rowMeans(ancom_decision_mat)
# ancom_roc_dec <- roc(difftaxa_truth, ancom_taxon_decision)
# ancom_difftaxa <- ancom_taxon_decision >= 0.9
# table(difftaxa_truth, ancom_difftaxa)

# ZOIB performance
zinf_decision <- combination_replicate$ZOIB_zinf_pval < 0.05
zinf_decision[is.na(zinf_decision)] <- FALSE
oinf_decision <- combination_replicate$ZOIB_oinf_pval < 0.05/2
oinf_decision[is.na(oinf_decision)] <- FALSE
beta_decision <- combination_replicate$ZOIB_beta_pval < 0.05/3
beta_decision[is.na(beta_decision)] <- FALSE

ZOIB_decision <- zinf_decision | oinf_decision | beta_decision
ZOIB_decision_mat <- map_significance(ZOIB_decision)
ZOIB_taxon_prob <- rowMeans(ZOIB_decision_mat)
# ZOIB_roc_dec <- roc(difftaxa_truth, ZOIB_taxon_decision)
# auc(ZOIB_roc_dec)


# GLMM performance
GLMM_decision <- combination_replicate$GLMM_pval < 0.05
GLMM_decision_mat <- map_significance(GLMM_decision)
GLMM_taxon_prob <- rowMeans(GLMM_decision_mat)


taxa_decisions <- cbind(ancom_taxon_prob,
                        ZOIB_taxon_prob,
                        GLMM_taxon_prob) %>% as.data.frame()
colnames(taxa_decisions) <- c("ANCOM", "ZOIB", "GLMM")
taxa_decisions$truth <- difftaxa_truth


threshold_selection <- function(rocobj){
  TP <- length(rocobj$cases) * rocobj$sensitivities
  FP <- length(rocobj$controls) * (1-rocobj$specificities)
  FDR <- FP / (FP + TP)
  cutoff <- rocobj$thresholds[min(which(FDR < 0.05))]
  return(cutoff)
}

ancom_roc_dec <- roc(taxa_decisions$truth, taxa_decisions$ANCOM)
ancom_threshold <- threshold_selection(ancom_roc_dec)
ancom_taxon_decision <- taxa_decisions$ANCOM > ancom_threshold
table(taxa_decisions$truth, ancom_taxon_decision)


ZOIB_roc_dec <- roc(taxa_decisions$truth, taxa_decisions$ZOIB)
ZOIB_threshold <- threshold_selection(ZOIB_roc_dec)
ZOIB_taxon_decision <- taxa_decisions$ZOIB > ZOIB_threshold
table(taxa_decisions$truth, ZOIB_taxon_decision)


GLMM_roc_dec <- roc(taxa_decisions$truth, taxa_decisions$GLMM)
GLMM_threshold <- threshold_selection(GLMM_roc_dec)
GLMM_taxon_decision <- taxa_decisions$GLMM > GLMM_threshold
table(taxa_decisions$truth, GLMM_taxon_decision)



# ancom_difftaxa <- ancom_taxon_decision >= 0.9


table(combination_replicate$diffratio, ancom_decision)
table(combination_replicate$diffratio, ZOIB_decision)
table(combination_replicate$diffratio, GLMM_decision)

# table(truth, GLMM_decision)
table(truth, ZOIB_decision)
# table(GLMM_decision, ZOIB_decision)
# table(ancom_decision, GLMM_decision)
table(ancom_decision, ZOIB_decision)

