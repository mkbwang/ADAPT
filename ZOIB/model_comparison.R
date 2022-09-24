rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
library(dplyr)


data <- readRDS(file.path(folder, 'simulated_data.rds'))
abn_info <- data$mean.eco.abn

taxa_names <- row.names(abn_info)
taxa_pairs <- combn(taxa_names, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("T1", "T2")

# collect the truth
truth <- taxa_pairs
truth$diffratio <- TRUE
for (j in 1:nrow(truth)){
  t1 <- truth$T1[j]
  t2 <- truth$T2[j]
  truth$diffratio[j] <- abn_info[t1, 'effect.size'] != abn_info[t2, 'effect.size']
}

# load model results
ANCOM_result <- readRDS(file.path(folder, 'ancom_result.rds'))
GLMM_result <- readRDS(file.path(folder, 'glmm_result.rds'))
ZOIB_result <- readRDS(file.path(folder, 'ZOIB_result.rds'))

combination <- cbind(truth, ANCOM_result$pval, GLMM_result$pval, ZOIB_result$zinfpval,
                     ZOIB_result$oinfpval, ZOIB_result$betapval)

colnames(combination) <- c("T1", "T2", "diffratio", "ANCOM_pval",
                           "GLMM_pval", 'ZOIB_zinf_pval', 'ZOIB_oinf_pval',
                           'ZOIB_beta_pval')


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


truth <- combination_replicate$diffratio
ancom_decision <- combination_replicate$ANCOM_pval <0.05
GLMM_decision <- combination_replicate$GLMM_pval <0.05

zinf_decision <- combination_replicate$ZOIB_zinf_pval < 0.05
zinf_decision[is.na(zinf_decision)] <- FALSE
oinf_decision <- combination_replicate$ZOIB_oinf_pval < 0.05
oinf_decision[is.na(oinf_decision)] <- FALSE
beta_decision <- combination_replicate$ZOIB_beta_pval < 0.05

ZOIB_decision <- zinf_decision | oinf_decision | beta_decision

table(truth, ancom_decision)
table(truth, GLMM_decision)
table(truth, beta_decision)


