rm(list=ls())

folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
library(dplyr)


# zoib_example <- read.csv(file.path(folder, 'ZOIB_result', 'zoib_result_70.csv'))

# fnum <- 1

# model_name <- 'glmdisp'

# Fisher's combination of p values, used for ZOIB
fisher_combine <- function(separatetests){
  pvalvec <- separatetests[!is.na(separatetests)]
  combined_statistic <- -2*sum(log(pvalvec))
  dof <- 2*length(pvalvec)
  combined_pval <- pchisq(combined_statistic, dof, lower.tail = FALSE)
  return(combined_pval)
}

# cutting threshold on the probability histogram
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

# map pairwise differential result onto a matrix
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


performance_summary <- function(fnum, model_name){
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

  file <- sprintf('%s_result_%d.csv', model_name, fnum)
  result <- read.csv(file.path(folder,
                                    sprintf('%s_result', model_name),
                                    file))

  if(model_name == 'zoib'){
    combination <- cbind(taxa_pairs, truth_vec,
                         result$zinfpval, result$oinfpval, result$betapval)
    colnames(combination) <- c("T1", "T2", "diffratio",
                               'zinfpval', 'oinfpval', 'betapval')
    combination$pval <- apply(combination[, c("zinfpval",
                                                                       "oinfpval",
                                                                       "betapval")],
                                             1, fisher_combine)
  } else{
    combination <- cbind(taxa_pairs, truth_vec, result$pval)
    colnames(combination) <- c("T1", "T2", "diffratio", 'pval')
  }

  combination_replicate <- combination
  combination_replicate$pval <- p.adjust(combination_replicate$pval,
                                              method='BH')

  pair_decision <- combination_replicate$pval < 0.05
  pair_decision_mat <- map_significance(pair_decision)
  taxon_prob <- rowMeans(pair_decision_mat)

  difftaxa_truth <- abn_info$effect.size != 1

  threshold <- threshold_selection2(taxon_prob)
  taxon_decision <- taxon_prob > threshold
  taxa_performance <- table(difftaxa_truth, taxon_decision) %>% as.vector()

  pair_performance <- table(combination$diffratio, pair_decision) %>% as.vector()


  overall_performance <- c(pair_performance, taxa_performance)
  return(list(taxon_prob = taxon_prob,
              threshold = threshold,
              performance = overall_performance))
}


glmdisp_result <- performance_summary(1, 'glmdisp')
ANCOM_result <- performance_summary(1, 'ancom')
zoib_result <- performance_summary(1, 'zoib')

taxa_probs <- cbind(ANCOM_result$taxon_prob,
                    zoib_result$taxon_prob,
                    glmdisp_result$taxon_prob) %>% as.data.frame()

glmdisp_threshold <- glmdisp_result$threshold
ANCOM_threshold <- ANCOM_result$threshold
zoib_threshold <- zoib_result$threshold

colnames(taxa_probs) <- c("ANCOM", "ZOIB", "GLMdisp")

library(ggplot2)


hist_ancom <- ggplot(taxa_probs, aes(x=ANCOM)) +
  geom_histogram(color="black", fill="white", binwidth = 0.05) +
  geom_vline(xintercept=ANCOM_threshold, linetype="dashed", color = "red")+
  xlab("W/(K-1)") + ylab("Taxa Count")+
  ggtitle("Wilkoxon")+xlim(0, 1)+
  theme_bw()

hist_ZOIB <- ggplot(taxa_probs, aes(x=ZOIB)) +
  geom_histogram(color="black", fill="white", binwidth = 0.05) +
  geom_vline(xintercept=zoib_threshold, linetype="dashed", color = "red")+
  xlab("W/(K-1)") + ylab("Taxa Count")+
  ggtitle("ZOIB")+xlim(0, 1)+
  theme_bw()

hist_GLMdisp <- ggplot(taxa_probs, aes(x=GLMdisp)) +
  geom_histogram(color="black", fill="white", binwidth = 0.05) +
  geom_vline(xintercept=glmdisp_threshold, linetype="dashed", color = "red")+
  xlab("W/(K-1)") + ylab("Taxa Count")+
  ggtitle("GLMdisp")+xlim(0, 1)+
  theme_bw()

library(cowplot)

hists <- plot_grid(hist_ancom,
                   hist_ZOIB,
                   hist_GLMdisp, ncol=3)

