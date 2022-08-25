# zero-one inflated beta regression

library(zoib)
library(dplyr)

#--------------following code was run on biostat cluster---------------#

# load all the data
rm(list=ls())
# data_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/RMD/CAARS_data'
# model_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/RMD/CAARS_Model_Summary'
folder = '/home/wangmk/MDAWG/DiffRatio'
counts_table <- read.csv(file.path(folder, 'filtered_count_order_level.csv'),
                         row.names = 1)
metadata <- read.csv(file.path(folder, 'filtered_metadata_order_level.csv'))

# identify pairs of taxa that generates different outcomes under different model setting
pvals_table <- read.csv(file.path(folder, 'Pvals_order_level.csv'))

# wilx_adjusted <- p.adjust(pvals_table$wilx)
# GLMM_adjusted <- p.adjust(pvals_table$GLMM_GQ_BOBYQA)

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

outcome <- foreach(j=1:nrow(pvals_table), .combine=rbind,
                       .packages=c('zoib', 'dplyr'),
                       .inorder=FALSE,
                        .errorhandling="remove") %dopar% {
     t1 <- pvals_table$Taxa1[j]
     t2 <- pvals_table$Taxa2[j]

     counts_t1 <- counts_table[, t1]
     counts_t2 <- counts_table[, t2]

     t1pt2n <- (counts_t1 > 0) & (counts_t2 == 0)
     t1nt2p <- (counts_t1 == 0) & (counts_t2 > 0)
     t1prop <- counts_t1 / (counts_t1 + counts_t2)

     df <- cbind(metadata$asthma, t1pt2n, t1nt2p, t1prop, counts_t1, counts_t2) %>% as.data.frame()
     colnames(df) <- c("Asthma", "T1PT2N", 'T1NT2P', 'T1prop', 'count_t1', 'count_t2')

     # first check if there is difference in zero inflation
     zinfmodel <- glm(T1NT2P ~ Asthma, data=df,
                      family = binomial(link = "logit"))
     zinfpval <- summary(zinfmodel)$coefficients[2,4]
     # pvals_table$infl0[j] <- zinfpval

     df_filtered <- df %>% filter(!T1NT2P)
     oinfmodel <- glm(T1PT2N ~ Asthma, data=df_filtered,
                      family = binomial(link = "logit"))
     oinfpval <- summary(oinfmodel)$coefficients[2,4]
     # pvals_table$infl1[j] <- oinfpval

     df_dp <- df %>% filter(count_t1 > 0 & count_t2 > 0)
     betapval <- NA
     if (nrow(df_dp) > 0){
       invisible(betamodel <- zoib(T1prop ~ Asthma|Asthma,
                         data=df_dp, random=0, zero.inflation = FALSE,
                         one.inflation = FALSE, joint=FALSE,
                         n.iter=10000, n.thin=20, n.burn=2000))
       betareg_summary <- summary(betamodel$coeff)$statistics
       beta_test_statistic <- -abs(betareg_summary[2,1] / betareg_summary[2,2])
       betapval <- pnorm(beta_test_statistic)
       # pvals_table$betareg[j] <- betapval
     }
     result = list(Rownum=j, Inf0 = zinfpval, Inf1=oinfpval,
                   betapval=betapval)
     return(as.data.frame(result))
}

write.csv(outcome, file.path(folder, 'Pvals_order_level_zoib.csv'),
          row.names = FALSE)
