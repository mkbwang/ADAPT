rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
# folder = '.'
library(dplyr)

fnum <- 1
input_fname <- sprintf('simulated_data_%d.rds', fnum)
data <- readRDS(file.path(folder, 'data', input_fname))

abn_info <- data$mean.eco.abn

taxa_names <- row.names(abn_info)
taxa_pairs <- combn(taxa_names, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("T1", "T2")
counts <- data$obs.abn

indicator <- data$grp - 1


# GLMM

library(lme4)
glmm_result <- taxa_pairs
glmm_result$effect <- 0
glmm_result$SE <- 0
glmm_result$pval <- 0

indicator <- data$grp - 1

library(foreach)
library(doParallel)
cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


ptm <- proc.time()
outcome <- foreach(j=1:nrow(glmm_result), .combine=rbind,
        .packages=c('dplyr', 'lme4'), .inorder=FALSE) %dopar%{
  result <- list(ID=j, effect=NA, SE=NA, pval=NA)
  t1 <- as.character(glmm_result$T1[j])
  t2 <- as.character(glmm_result$T2[j])
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


glmm_result$effect[outcome$ID] <- outcome$effect
glmm_result$SE[outcome$ID] <- outcome$SE
glmm_result$pval[outcome$ID] <- outcome$pval


filename <- sprintf('glmm_result_%d.csv', fnum)
write.csv(ancom_result, file.path(folder, 'glmm_result', filename),
          row.names = FALSE)


