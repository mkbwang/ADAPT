rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
library(dplyr)

data <- readRDS(file.path(folder, 'simulated_data.rds'))
abn_info <- data$mean.eco.abn
indicator <- data$grp - 1

taxa_names <- row.names(abn_info)
taxa_pairs <- combn(taxa_names, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("T1", "T2")

counts <- data$obs.abn


# GLMM

library(lme4)
glmm_result <- taxa_pairs
glmm_result$effect <- 0
glmm_result$SE <- 0
glmm_result$pval <- 0

ptm <- proc.time()
outcome <- foreach(j=1:nrow(glmm_result), .combine=rbind,
        .packages=c('dplyr', 'lme4'), .inorder=FALSE,
        .errorhandling="remove") %dopar%{
  result <- list(ID=j, effect=NA, SE=NA, pval=NA)
  t1 <- glmm_result$T1[j]
  t2 <- glmm_result$T2[j]
  selected_counts <- counts[c(t1, t2),] %>%
    t() %>% as.data.frame()
  colnames(selected_counts) <- c('T1', 'T2')
  selected_counts$covar <- indicator
  selected_counts$ID <- row.names(selected_counts)
  mixed_model <- suppressMessages(glmer(cbind(T1,T2) ~ covar+(1|ID), family="binomial", data=selected_counts,
                                        nAGQ=20)) %>%
    summary()
  result$effect <- mixed_model$coefficients[2, 1]
  result$SE <- mixed_model$coefficients[2, 2]
  result$pval <- mixed_model$coefficients[2, 4]
  return(as.data.frame(result))
}
duration <- proc.time() - ptm
glmm_result$effect <- outcome$effect
glmm_result$SE <- outcome$SE
glmm_result$pval <- outcome$pval

saveRDS(glmm_result, file=file.path(folder, 'glmm_result.rds'))


