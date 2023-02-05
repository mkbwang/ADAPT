rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'

library(dplyr)
fnum <- 1
input_fname <- sprintf('simulated_data_%d.rds', fnum)
data <- readRDS(file.path(folder, 'data', 'semiparametric', input_fname))

counts <- data$otu.tab.sim

taxa_names <- row.names(counts)
taxa_pairs <- combn(taxa_names, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("T1", "T2")

indicator <- as.vector(data$covariate)

# try ANCOM
ancom_result <- taxa_pairs
ancom_result$pval <- NA


library(foreach)
library(doParallel)
cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


ptm <- proc.time()
outcome <- foreach(j=1:nrow(ancom_result), .combine=rbind,
                   .packages=c('dplyr'), .inorder=FALSE,
                   .errorhandling="remove") %dopar% {
  t1 <- as.character(ancom_result$T1[j])
  t2 <- as.character(ancom_result$T2[j])
  selected_counts <- counts[c(t1, t2),] %>%
    t() %>% as.data.frame()
  colnames(selected_counts) <- c('T1', 'T2')
  selected_counts$covar <- indicator
  selected_counts$logratio <- log(selected_counts$T1+1) - log(selected_counts$T2+1)
  logratio_test <- wilcox.test(logratio~covar, data=selected_counts)
  return(c(j, logratio_test$p.value))
}
duration <- proc.time() - ptm



ancom_result$pval[outcome[, 1]] <- outcome[, 2]

filename <- sprintf('ancom_result_%d.csv', fnum)
write.csv(ancom_result, file.path(folder, 'result', filename),
          row.names = FALSE)


