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

GLMM_result <- readRDS(file.path(folder, 'glmm_result.rds'))

combination <- cbind(truth, ANCOM_result$pval, ZOIB_result$zinfpval,
                     ZOIB_result$oinfpval, ZOIB_result$betapval)
colnames(combination) <- c("T1", "T2", "diffratio", "ANCOM_pval",
                            'ZOIB_zinf_pval', 'ZOIB_oinf_pval',
                           'ZOIB_beta_pval')


combination_replicate <- combination

combination_replicate$ANCOM_pval <- p.adjust(combination_replicate$ANCOM_pval,
                                             method='BH')
# combination_replicate$GLMM_pval <- p.adjust(combination_replicate$GLMM_pval,
#                                             method='BH')
combination_replicate$ZOIB_zinf_pval <- p.adjust(combination_replicate$ZOIB_zinf_pval,
                                            method='BH')
combination_replicate$ZOIB_oinf_pval <- p.adjust(combination_replicate$ZOIB_oinf_pval,
                                                 method='BH')
combination_replicate$ZOIB_beta_pval <- p.adjust(combination_replicate$ZOIB_beta_pval,
                                                 method='BH')


truth <- combination_replicate$diffratio
ancom_decision <- combination_replicate$ANCOM_pval <0.05

map_significance <- function(pairresult){
  decision_mat <- matrix(FALSE, nrow=1000, ncol=1000)
  cursor <- 1
  for (j in 1:999){
    decision_mat[j, (j+1):1000] <- pairresult[cursor: (cursor+1000-j-1)]
    cursor <- cursor + 1000-j
  }
  decision_mat <- decision_mat + t(decision_mat)
  return(decision_mat)
}

# GLMM_decision <- combination_replicate$GLMM_pval <0.05
difftaxa_truth <- abn_info$effect.size != 1

library(pROC)
ancom_decision_mat <- map_significance(ancom_decision)
ancom_taxon_decision <- rowMeans(ancom_decision_mat)
ancom_roc_dec <- roc(difftaxa_truth, ancom_taxon_decision)
auc(ancom_roc_dec)
# ancom_difftaxa <- ancom_taxon_decision >= 0.9
# table(difftaxa_truth, ancom_difftaxa)


zinf_decision <- combination_replicate$ZOIB_zinf_pval < 0.05
zinf_decision[is.na(zinf_decision)] <- FALSE
oinf_decision <- combination_replicate$ZOIB_oinf_pval < 0.05
oinf_decision[is.na(oinf_decision)] <- FALSE
beta_decision <- combination_replicate$ZOIB_beta_pval < 0.05
beta_decision[is.na(beta_decision)] <- FALSE

ZOIB_decision <- zinf_decision | oinf_decision | beta_decision
ZOIB_decision_mat <- map_significance(ZOIB_decision)
ZOIB_taxon_decision <- rowMeans(ZOIB_decision_mat)
ZOIB_roc_dec <- roc(difftaxa_truth, ZOIB_taxon_decision)
auc(ZOIB_roc_dec)


table(truth, ancom_decision)
# table(truth, GLMM_decision)
table(truth, ZOIB_decision)
# table(GLMM_decision, ZOIB_decision)
# table(ancom_decision, GLMM_decision)
table(ancom_decision, ZOIB_decision)


combination$ANCOM_decision <- ancom_decision
combination$GLMM_decision <- GLMM_decision
combination$ZOIB_decision <- ZOIB_decision


pairs1 <- combination %>% filter(!diffratio & !GLMM_decision & ZOIB_decision)
pairs2 <- combination %>% filter(diffratio & !GLMM_decision & ZOIB_decision)

sample_mat <- data$obs.abn


## example 1



