# zero-one inflated beta regression

library(zoib)
library(dplyr)

#--------------following code was run on biostat cluster---------------#

# load all the data
rm(list=ls())
data_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/CAARS/CAARS_data'
model_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/CAARS/CAARS_Model_Summary'
# data_folder = '/home/wangmk/MDAWG/DiffRatio'
# model_folder = '/home/wangmk/MDAWG/DiffRatio'
counts_table <- read.csv(file.path(data_folder, 'filtered_count_order_level.csv'),
                         row.names = 1)
metadata <- read.csv(file.path(data_folder, 'filtered_metadata_order_level.csv'))

# identify pairs of taxa that generates different outcomes under different model setting
pvals_table <- read.csv(file.path(model_folder, 'Pvals_order_level.csv'))

j <- 87

t1 <- pvals_table$Taxa1[j]
t2 <- pvals_table$Taxa2[j]

counts_t1 <- counts_table[, t1]
counts_t2 <- counts_table[, t2]

t1pt2n <- (counts_t1 > 0) & (counts_t2 == 0)
t1nt2p <- (counts_t1 == 0) & (counts_t2 > 0)
t1nt2n <- (counts_t1 == 0) & (counts_t2 == 0)
t1prop <- counts_t1 / (counts_t1 + counts_t2)

df <- cbind(metadata$asthma, t1pt2n, t1nt2p, t1nt2n, t1prop, counts_t1, counts_t2) %>% as.data.frame()
colnames(df) <- c("Asthma", "T1PT2N", 'T1NT2P', 'T1NT2N', 'T1prop', 'count_t1', 'count_t2')
df_permuted <- df
df_permuted$Asthma <- sample(df_permuted$Asthma)

df <- df %>% filter(T1NT2N == 0)
df_permuted <- df_permuted %>% filter(T1NT2N == 0)

ZOIB_test <- function(mydf){
        # first check if there is difference in zero inflation
        zinfmodel <- glm(T1NT2P ~ Asthma, data=mydf,
                         family = binomial(link = "logit"))
        zinfpval <- summary(zinfmodel)$coefficients[2,4]
        # pvals_table$infl0[j] <- zinfpval

        mydf_filtered <- mydf %>% filter(!T1NT2P)
        oinfmodel <- glm(T1PT2N ~ Asthma, data=mydf_filtered,
                         family = binomial(link = "logit"))
        oinfpval <- summary(oinfmodel)$coefficients[2,4]
        # pvals_table$infl1[j] <- oinfpval

        mydf_dp <- mydf %>% filter(count_t1 > 0 & count_t2 > 0)
        betapval <- NA
        if (nrow(mydf_dp) > 0){
                invisible(betamodel <- zoib(T1prop ~ Asthma|Asthma,
                                            data=mydf_dp, random=0, zero.inflation = FALSE,
                                            one.inflation = FALSE, joint=FALSE,
                                            n.iter=10000, n.thin=20, n.burn=2000))
                betareg_summary <- summary(betamodel$coeff)$statistics
                beta_test_statistic <- -abs(betareg_summary[2,1] / betareg_summary[2,2])
                betapval <- pnorm(beta_test_statistic)
                # pvals_table$betareg[j] <- betapval
        }
        result = list(Inf0 = zinfpval, Inf1=oinfpval,
                      betapval=betapval)
        return(result)
}


original_testresult <- ZOIB_test(df)
permuted_testresult <- ZOIB_test(df_permuted)

originaldata_result <- sprintf("%s, %s, %f, %f, %f", t1, t2, original_testresult$Inf0,
                               original_testresult$Inf1, original_testresult$betapval)
permuteddata_result <- sprintf("%s, %s, %f, %f, %f", t1, t2, permuted_testresult$Inf0,
                               permuted_testresult$Inf1, permuted_testresult$betapval)

write(originaldata_result, file=file.path(model_folder, 'Pvals_order_level_ZOIB.csv'),
      append=TRUE)
write(permuteddata_result, file=file.path(model_folder, 'Pvals_order_level_ZOIB_permuted.csv'),
      append=TRUE)
