# Check out the p value and effect size distributions for all the pairs of taxa
library(dplyr)
data_folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data'

pairGLM_folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/glmdisp_result'

replicate = 1

# load simulated data
simulated_data = readRDS(file.path(data_folder,
                                   sprintf("simulated_data_%d.rds", replicate)))

true_abs_abundance <- simulated_data$mean.eco.abn
true_abs_abundance$effect.size <- true_abs_abundance$temp.grp2 / true_abs_abundance$temp.grp1

# high_increase_taxa <- true_abs_abundance %>% filter(effect.size > 2)
# low_increase_taxa <- true_abs_abundance %>% filter(effect.size > 1 & effect.size < 2)
# stable_taxa <- true_abs_abundance %>% filter(effect.size == 1)
# low_decrease_taxa <- true_abs_abundance %>% filter(effect.size < 1 & effect.size > 0.5)
# high_decrease_taxa <- true_abs_abundance %>% filter(effect.size <0.5)

set.seed(2020)
# example_taxon_highincrease <- sample(row.names(high_increase_taxa), size=1)
# high_increase_taxa[example_taxon_highincrease, ]
#
# example_taxon_lowincrease <- sample(row.names(low_increase_taxa), size=1)
# low_increase_taxa[example_taxon_lowincrease, ]
#
# example_taxon_stable <- sample(row.names(stable_taxa), size=1)
# stable_taxa[example_taxon_stable, ]
#
# example_taxon_lowdecrease <- sample(row.names(low_decrease_taxa), size=1)
# low_decrease_taxa[example_taxon_lowdecrease, ]
#
# example_taxon_highdecrease <- sample(row.names(high_decrease_taxa), size=1)
# high_decrease_taxa[example_taxon_highdecrease, ]


# load GLM result
GLM_result = read.csv(file.path(pairGLM_folder,
                                sprintf("glmdisp_result_%d.csv", replicate)))


# GLM_result_highincrease <- GLM_result %>% filter(T1 == example_taxon_highincrease |
#                                                    T2 == example_taxon_highincrease)
# GLM_result_lowincrease <-  GLM_result %>% filter(T1 == example_taxon_lowincrease |
#                                                    T2 == example_taxon_lowincrease)
GLM_result_stable <- GLM_result %>% filter(T1 == "taxon2" |
                                             T2 == "taxon2")

# GLM_result_lowdecrease <-  GLM_result %>% filter(T1 == example_taxon_lowdecrease |
#                                                    T2 == example_taxon_lowdecrease)
# GLM_result_highdecrease <-  GLM_result %>% filter(T1 == example_taxon_highdecrease |
#                                                    T2 == example_taxon_highdecrease)

uniq_taxon <- sprintf("taxon%d", seq(1, 1000))

# remove extreme outliers to estimate normal distribution parameter correctly
IOR <- function(myvec){
  count <- 0
  finished <- FALSE
  while(!finished & count < 20){
    avg <- mean(myvec)
    sdv <- sd(myvec)
    upperbound <- avg + 3*sdv
    lowerbound <- avg - 3*sdv
    newvec <- myvec[myvec < upperbound & myvec > lowerbound]
    if(length(newvec) == length(myvec)){
      myvec <- newvec
      finished=TRUE
    } else{
      myvec <- newvec
      count <- count + 1
    }
  }
  return(myvec)
}

# test statistics regarding a specific taxon and its normal distribution estimate
extract_teststat <- function(glmresult, taxon){
  fold_change <- true_abs_abundance[taxon, 'effect.size']
  glm_result_subset <- glmresult %>% filter(T1 == taxon |
                                                 T2 == taxon)
  teststat <- glm_result_subset$effect / glm_result_subset$SE
  teststat[glm_result_subset$T2 == taxon] <- -teststat[glm_result_subset$T2 == taxon]
  filtered_teststat <- IOR(teststat)
  fitnormal <- fitdistr(filtered_teststat, densfun="normal")
  return(list(teststat=filtered_teststat, estimate=fitnormal$estimate))
}

# example <- extract_teststat(GLM_result, "taxon437")

library(foreach)
library(doParallel)
cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


summary_teststat <- foreach(i=1:length(uniq_taxon), .combine=rbind,
                            .packages=c("dplyr", "MASS"), .inorder=FALSE) %dopar%{
  taxon_characteristics <- extract_teststat(GLM_result, taxon=uniq_taxon[i])
  return(taxon_characteristics$estimate)
}


summary_teststat <- as.data.frame(summary_teststat)
colnames(summary_teststat) <- c("Mean", "SD")
row.names(summary_teststat) <- sprintf("Taxon%d", seq(1, 1000))

library(philentropy)
# randommat <- matrix(runif(8, 0, 8), nrow=4, ncol=2)

hellinger_dist_norm <- function(params1, params2){
  mu1 <- params1[1]
  sigma1 <- params1[2]
  mu2 <- params2[1]
  sigma2 <- params2[2]
  sqHdist <- 1-sqrt(2*sigma1*sigma2/(sigma1^2+sigma2^2))*
    exp(-0.25*(mu1-mu2)^2/(sigma1^2+sigma2^2))

  return(sqrt(sqHdist))
}

distmat <- proxy::dist(summary_teststat, method=hellinger_dist_norm)

hier_clust <- hclust(distmat, method="average")

output <- cutree(hclust(distmat), k=3)
pred_label <- c("Stable", "Increase", "Decrease")

true_abs_abundance$prediction <- output
true_abs_abundance$category <- (true_abs_abundance$effect.size == 1)+
  3*(true_abs_abundance$effect.size < 1)+
  2*(true_abs_abundance$effect.size > 1)

table(true_abs_abundance$prediction, true_abs_abundance$category)


# hclusCut <- function(x, k, ...){
#   distmat <- proxy::dist(x, method=hellinger_dist_norm)
#   list(cluster = cutree(hclust(distmat), k=k))
# }
#
# hier_clust <- clusGap(summary_teststat, FUNcluster = hclusCut, K.max=15, B=60)

# stat_highincrease <- extract_teststat(GLM_result_highincrease, example_taxon_highincrease)
# stat_lowincrease <- extract_teststat(GLM_result_lowincrease, example_taxon_lowincrease)
# stat_stable <- extract_teststat(GLM_result_stable, example_taxon_stable)
# stat_lowdecrease <- extract_teststat(GLM_result_lowdecrease, example_taxon_lowdecrease)
# stat_highdecrease <- extract_teststat(GLM_result_highdecrease, example_taxon_highdecrease)
#
#
# stat_combined <- rbind(stat_highincrease, stat_lowincrease,
#                        stat_stable, stat_lowdecrease,
#                        stat_highdecrease)
#
# library(ggplot2)
#
# ggplot(stat_combined, aes(x=teststat)) +
#   geom_histogram(color="black", fill="white", binwidth=0.5) +
#   facet_wrap(~taxon, ncol=1) + xlim(-20, 20)+ geom_vline(xintercept=0, linetype=2)+
#   xlab("Test Statistic") + ylab("Count") + theme_bw()


library(MASS)

# highincrease_fit <- fitdistr(stat_highincrease$teststat, "t",
#          start = list(m = median(stat_highincrease$teststat),
#                       s = sd(stat_highincrease$teststat), df=2), lower=c(0.001, 0.001, 0.001))


# highincrease_fit <- fitdistr(stat_highincrease$teststat, 'normal')
# lowincrease_fit <- fitdistr(stat_lowincrease$teststat, 'normal')
# stable_fit <- fitdistr(stat_stable$teststat, 'normal')
# lowdecrease_fit <- fitdistr(stat_lowdecrease$teststat, 'normal')
# highdecrease_fit <- fitdistr(stat_highdecrease$teststat, 'normal')



diff1 <- hellinger_dist_norm(highdecrease_fit, lowdecrease_fit)


