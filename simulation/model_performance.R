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

fisher_combine <- function(separatetests){
  pvalvec <- separatetests[!is.na(separatetests)]
  combined_statistic <- -2*sum(log(pvalvec))
  dof <- 2*length(pvalvec)
  combined_pval <- pchisq(combined_statistic, dof, lower.tail = FALSE)
  return(combined_pval)
}

combination_replicate$ZOIB_pval <- apply(combination_replicate[, c("ZOIB_zinf_pval",
                                                                   "ZOIB_oinf_pval",
                                                                   "ZOIB_beta_pval")],
                                         1, fisher_combine)


combination_replicate$ZOIB_pval <- p.adjust(combination_replicate$ZOIB_pval,
                                            method="BH")


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

# library(pROC)
ancom_decision_mat <- map_significance(ancom_decision)
ancom_taxon_prob <- rowMeans(ancom_decision_mat)


## ZOIB performance

ZOIB_decision <- combination_replicate$ZOIB_pval < 0.05
ZOIB_decision_mat <- map_significance(ZOIB_decision)
ZOIB_taxon_prob <- rowMeans(ZOIB_decision_mat)



## GLMM performance
GLMM_decision <- combination_replicate$GLMM_pval < 0.05
GLMM_decision_mat <- map_significance(GLMM_decision)
GLMM_taxon_prob <- rowMeans(GLMM_decision_mat)


taxa_decisions <- cbind(ancom_taxon_prob,
                        ZOIB_taxon_prob,
                        GLMM_taxon_prob) %>% as.data.frame()
colnames(taxa_decisions) <- c("ANCOM", "ZOIB", "GLMM")
taxa_decisions$truth <- difftaxa_truth


threshold_selection1 <- function(rocobj){
  TP <- length(rocobj$cases) * rocobj$sensitivities
  FP <- length(rocobj$controls) * (1-rocobj$specificities)
  FDR <- FP / (FP + TP)
  cutoff <- rocobj$thresholds[min(which(FDR < 0.05))]
  return(cutoff)
}

threshold_selection2 <- function(probvecs){
  pop_mean <- mean(probvecs)
  uniqvals <- unique(probvecs) %>% sort()
  j <- 1
  group1 <- probvecs[probvecs <= uniqvals[j]]
  group1mean <- mean(group1)
  group2 <- probvecs[probvecs > uniqvals[j]]
  group2mean <- mean(group2)
  # seek largest intercluster variance
  intercluster <- length(group1) * (group1mean - pop_mean)^2 +
    length(group2) * (group2mean - pop_mean)^2
  while(1){
    group1 <- probvecs[probvecs <= uniqvals[j+1]]
    group1mean <- mean(group1)
    group2 <- probvecs[probvecs > uniqvals[j+1]]
    group2mean <- mean(group2)
    newintercluster <- length(group1) * (group1mean - pop_mean)^2 +
      length(group2) * (group2mean - pop_mean)^2
    if (newintercluster > intercluster){
      intercluster <- newintercluster
      j <- j+1
    } else{
      break
    }
  }
  return(uniqvals[j])
}



# ancom_roc_dec <- roc(taxa_decisions$truth, taxa_decisions$ANCOM)
ancom_threshold <- threshold_selection2(ancom_taxon_prob)
ancom_taxon_decision <- taxa_decisions$ANCOM > ancom_threshold
taxa_ancom_performance <- table(taxa_decisions$truth, ancom_taxon_decision) %>% as.vector()

library(ggplot2)


# ZOIB_roc_dec <- roc(taxa_decisions$truth, taxa_decisions$ZOIB)
ZOIB_threshold <- threshold_selection2(ZOIB_taxon_prob)
ZOIB_taxon_decision <- taxa_decisions$ZOIB > ZOIB_threshold
taxa_ZOIB_performance <- table(taxa_decisions$truth, ZOIB_taxon_decision) %>% as.vector()


# GLMM_roc_dec <- roc(taxa_decisions$truth, taxa_decisions$GLMM)
GLMM_threshold <- threshold_selection2(GLMM_taxon_prob)
GLMM_taxon_decision <- taxa_decisions$GLMM > GLMM_threshold
taxa_GLMM_performance <- table(taxa_decisions$truth, GLMM_taxon_decision) %>% as.vector()


hist_ancom <- ggplot(taxa_decisions, aes(x=ANCOM)) +
  geom_histogram(color="black", fill="white", binwidth = 0.05) +
  geom_vline(xintercept=ancom_threshold, linetype="dashed", color = "red")+
  xlab("W/(K-1)") + ylab("Taxa Count")+
  ggtitle("Wilkoxon")+xlim(0, 1)+
  theme_bw()

hist_ZOIB <- ggplot(taxa_decisions, aes(x=ZOIB)) +
  geom_histogram(color="black", fill="white", binwidth = 0.05) +
  geom_vline(xintercept=ZOIB_threshold, linetype="dashed", color = "red")+
  xlab("W/(K-1)") + ylab("Taxa Count")+
  ggtitle("ZOIB")+xlim(0, 1)+
  theme_bw()

hist_GLMM <- ggplot(taxa_decisions, aes(x=GLMM)) +
  geom_histogram(color="black", fill="white", binwidth = 0.05) +
  geom_vline(xintercept=GLMM_threshold, linetype="dashed", color = "red")+
  xlab("W/(K-1)") + ylab("Taxa Count")+
  ggtitle("GLMM")+xlim(0, 1)+
  theme_bw()

library(cowplot)

hists <- plot_grid(hist_ancom,
                   hist_ZOIB,
                   hist_GLMM, ncol=3)

pair_ancom_performance <- table(combination_replicate$diffratio, ancom_decision) %>% as.vector()
pair_ZOIB_performance <- table(combination_replicate$diffratio, ZOIB_decision) %>% as.vector()
pair_GLMM_performance <- table(combination_replicate$diffratio, GLMM_decision) %>% as.vector()


ancom_performance <- c(pair_ancom_performance, taxa_ancom_performance)
ZOIB_performance <- c(pair_ZOIB_performance, taxa_ZOIB_performance)
GLMM_performance <- c(pair_GLMM_performance, taxa_GLMM_performance)

write(paste(ancom_performance, collapse=','),
      file=file.path(folder, 'ancom_result', 'ancom_performance.csv'),
      append=TRUE)
write(paste(ZOIB_performance, collapse=','),
      file=file.path(folder, 'ZOIB_result', 'ZOIB_performance.csv'),
      append=TRUE)
write(paste(GLMM_performance, collapse=','),
      file=file.path(folder, 'glmm_result', 'glmm_performance.csv'),
      append=TRUE)
