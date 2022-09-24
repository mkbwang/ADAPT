# zero-one inflated beta regression
library(betareg)
library(dplyr)

#--------------following code was run on biostat cluster---------------#

# load all the data
rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'

data <- readRDS(file.path(folder, 'simulated_data.rds'))
abn_info <- data$mean.eco.abn
indicator <- data$grp - 1

taxa_names <- row.names(abn_info)
taxa_pairs <- combn(taxa_names, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("T1", "T2")

# start zero.one inflated beta regression
zoib_result <- taxa_pairs
zoib_result$zinfeffect <- NA
zoib_result$zinfpval <- NA
zoib_result$oinfeffect <- NA
zoib_result$oinfpval <- NA
zoib_result$betaeffect <- NA
zoib_result$betapval <- NA

counts <- data$obs.abn

ZOIB_test <- function(mydf){
  mydf <- mydf %>% filter(!t1nt2n)
  result = list(zinfeffect=NA, zinfpval=NA,
                oinfeffect=NA, oinfpval=NA,
                betaeffect=NA, betapval=NA)
  # first check if there is difference in zero inflation
  if (sum(mydf$t1nt2p) > 0){
    zinfmodel <- glm(t1nt2p ~ covar, data=mydf,
                     family = binomial(link = "logit")) %>% summary()
    result$zinfeffect <- zinfmodel$coefficients[2, 1]
    result$zinfpval <- zinfmodel$coefficients[2,4]
  }

  mydf_filtered <- mydf %>% filter(!t1nt2p)
  if (sum(mydf_filtered$t1pt2n) > 0){
    oinfmodel <- glm(t1pt2n ~ covar, data=mydf_filtered,
                     family = binomial(link = "logit")) %>% summary()
    result$oinfeffect <- oinfmodel$coefficients[2,1]
    result$oinfpval <- oinfmodel$coefficients[2,4]
  }

  mydf_dp <- mydf %>% filter(T1 > 0 & T2 > 0)
  if (nrow(mydf_dp) > 2){
    invisible(beta_model <- betareg(prop ~ covar | covar, data=mydf_dp) %>% summary())
    result$betaeffect <- beta_model$coefficients$mean[2, 1]
    result$betapval <- beta_model$coefficients$mean[2, 4]
  }

  return(result)
}


library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


ptm <- proc.time()
outcome <- foreach(j=1:nrow(zoib_result), .combine=rbind,
                   .packages=c('dplyr', 'betareg'), .inorder=FALSE) %dopar% {

   t1 <- zoib_result$T1[j]
   t2 <- zoib_result$T2[j]

   selected_counts <- counts[c(t1, t2), ] %>% t() %>%
     as.data.frame()
   colnames(selected_counts) <- c('T1', 'T2')
   selected_counts$covar <- indicator

   selected_counts$t1pt2n <- (selected_counts$T1 > 0) & (selected_counts$T2 == 0)
   selected_counts$t1nt2p <- (selected_counts$T1 == 0) & (selected_counts$T2 > 0)
   selected_counts$t1nt2n <- (selected_counts$T1 == 0) & (selected_counts$T2 == 0)
   selected_counts$prop <- selected_counts$T1/(selected_counts$T1 + selected_counts$T2)
   test_result <- ZOIB_test(selected_counts)
   test_result$ID <- j
  return(as.data.frame(test_result))
}
duration <- proc.time() - ptm

zoib_result$zinfeffect <- outcome$zinfeffect
zoib_result$zinfpval <- outcome$zinfpval
zoib_result$oinfeffect <- outcome$oinfeffect
zoib_result$oinfpval <- outcome$oinfpval
zoib_result$betaeffect <- outcome$betaeffect
zoib_result$betapval <- outcome$betapval

saveRDS(zoib_result, file=file.path(folder, 'ZOIB_result.rds'))

