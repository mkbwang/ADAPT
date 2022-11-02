# check out the FN and FP of each pairwise method
rm(list=ls())

folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
library(dplyr)


fnum <- 81
data_name <- sprintf('simulated_data_%d.rds', fnum)
data <- readRDS(file.path(folder, 'data', data_name))
abn_info <- data$mean.eco.abn
counts <- data$obs.abn

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


ancom_decision <- combination_replicate$ANCOM_pval <0.05
ZOIB_decision <- combination_replicate$ZOIB_pval < 0.05
GLMM_decision <- combination_replicate$GLMM_pval < 0.05

pairs_conclusion <- data.frame(ID = 1:length(ancom_decision),
  truth = truth_vec,
                               ancom = ancom_decision,
                               zoib = ZOIB_decision,
                               glmm = GLMM_decision)

library(ggplot2)
library("scales")

biplot_viz <- function(data, tpair){

  data$Proportion <- data[, tpair[1]] / (data[, tpair[1]] + data[, tpair[2]])
  data$logratio <- log(data[, tpair[1]] + 1) - log(data[, tpair[2]] + 1)

  count_plot <- ggplot(data, aes_string(x=tpair[2], y=tpair[1], color='Group')) +
    geom_jitter(alpha=0.6, width = 0.1, height = 0.1) + scale_x_continuous(breaks = pretty_breaks())+
    scale_y_continuous(breaks = pretty_breaks())+
    scale_color_manual(values=c("#E69F00", "#56B4E9"))+
    theme_bw()

  proportion_plot <- ggplot(data, aes(x=Proportion, color=Group, fill=Group)) +
    geom_histogram(position="identity", alpha=0.5)+
    scale_x_continuous(breaks = pretty_breaks())+
    scale_y_continuous(breaks = pretty_breaks())+
    scale_color_manual(values=c("#E69F00", "#56B4E9"))+
    xlim(-0.01, 1.01)+ scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    xlab(sprintf('%s/(%s+%s)', tpair[1], tpair[1], tpair[2])) + ylab("Count") +  theme_bw()

  logratio_plot <- ggplot(data, aes(x=Group, y=logratio)) +
    geom_violin() + geom_boxplot(width=0.1) +
    ylab(sprintf('log(%s+1)-log(%s+1)', tpair[1], tpair[2])) + xlab("Group")+
    coord_flip() + theme_bw()

  return(list(countplot = count_plot,
              propplot = proportion_plot,
              logratioplot = logratio_plot))

}


# wilcoxon FP

wilcox_FP <- pairs_conclusion %>% filter(!truth & ancom & !zoib & !glmm)
example_pair1 <- c(taxa_pairs$T1[wilcox_FP$ID[10]], taxa_pairs$T2[wilcox_FP$ID[10]])
example_counts1 <- counts[example_pair1, ] %>% t() %>% as.data.frame()
# example_counts1$proportion <- example_counts1$taxon1/(example_counts1$taxon1 + example_counts1$taxon91)
example_counts1$Group <- rep(c(0, 1), each=50)
example_counts1$Group <- as.factor(example_counts1$Group)

plot_wilcox_FP <- biplot_viz(example_counts1, example_pair1)
library(cowplot)
plot_wilcox_FP_combined <- plot_grid(plot_wilcox_FP$countplot,
          plot_wilcox_FP$propplot,
          plot_wilcox_FP$logratioplot,
          ncol=2)

abn_info[example_pair1, ]


# wilcoxon FN
wilcox_FN <- pairs_conclusion %>% filter(truth & !ancom & zoib & glmm)
example_pair2 <- c(taxa_pairs$T1[wilcox_FN$ID[8]],
                   taxa_pairs$T2[wilcox_FN$ID[8]])
example_counts2 <- counts[example_pair2, ] %>% t() %>% as.data.frame()
example_counts2$Group <- rep(c(0, 1), each=50)
example_counts2$Group <- as.factor(example_counts1$Group)
plot_wilcox_FN <- biplot_viz(example_counts2, example_pair2)

plot_wilcox_FN_combined <- plot_grid(plot_wilcox_FN$countplot,
          plot_wilcox_FN$propplot,
          plot_wilcox_FN$logratioplot,
          ncol=2)

abn_info[example_pair2, ]

ancom_feedback <- ANCOM_result[wilcox_FN$ID, ]
ZOIB_feedback <- ZOIB_result[wilcox_FN$ID, ]
GLMM_feedback <- GLMM_result[wilcox_FN$ID, ]

# ZOIB FP
ZOIB_FP <- pairs_conclusion %>% filter(!truth & !ancom & zoib & !glmm)
example_pair3 <- c(taxa_pairs$T1[ZOIB_FP$ID[5]],
                   taxa_pairs$T2[ZOIB_FP$ID[5]])
example_counts3 <- counts[example_pair3, ] %>% t() %>% as.data.frame()
example_counts3$Group <- rep(c(0, 1), each=50)
example_counts3$Group <- as.factor(example_counts1$Group)

plot_ZOIB_FP<- biplot_viz(example_counts3, example_pair3)

plot_ZOIB_FP_combined <- plot_grid(plot_ZOIB_FP$countplot,
          plot_ZOIB_FP$propplot,
          plot_ZOIB_FP$logratioplot,
          ncol=2)

abn_info[example_pair3, ]

ancom_feedback <- ANCOM_result[ZOIB_FP$ID[5], ]
ZOIB_feedback <- ZOIB_result[ZOIB_FP$ID[5], ]
GLMM_feedback <- GLMM_result[ZOIB_FP$ID[5], ]


# ZOIB FN
ZOIB_FN <- pairs_conclusion %>% filter(truth & ancom & !zoib & glmm)
example_pair4 <- c(taxa_pairs$T1[3779],
                   taxa_pairs$T2[3779])
example_counts4 <- counts[example_pair4, ] %>% t() %>% as.data.frame()
example_counts4$Group <- rep(c(0, 1), each=50)
example_counts4$Group <- as.factor(example_counts1$Group)

plot_ZOIB_FN <- biplot_viz(example_counts4, example_pair4)

plot_ZOIB_FN_combined <- plot_grid(plot_ZOIB_FN$countplot,
          plot_ZOIB_FN$propplot,
          plot_ZOIB_FN$logratioplot,
          ncol=2)

abn_info[example_pair4, ]

ancom_feedback <- ANCOM_result[3779, ]
ZOIB_feedback <- ZOIB_result[3779, ]
GLMM_feedback <- GLMM_result[3779, ]


# GLMM FP
GLMM_FP <- pairs_conclusion %>% filter(!truth & !ancom & !zoib & glmm)
example_pair5 <- c(taxa_pairs$T1[GLMM_FP$ID[6]],
                   taxa_pairs$T2[GLMM_FP$ID[6]])
example_counts5 <- counts[example_pair5, ] %>% t() %>% as.data.frame()
example_counts5$Group <- rep(c(0, 1), each=50)
example_counts5$Group <- as.factor(example_counts5$Group)


plot_GLMM_FP <- biplot_viz(example_counts5, example_pair5)
plot_GLMM_FP_combined <- plot_grid(plot_GLMM_FP$countplot,
          plot_GLMM_FP$propplot,
          plot_GLMM_FP$logratioplot,
          ncol=2)

ancom_feedback <- ANCOM_result[GLMM_FP$ID[6], ]
ZOIB_feedback <- ZOIB_result[GLMM_FP$ID[6], ]
GLMM_feedback <- GLMM_result[GLMM_FP$ID[6], ]

abn_info[example_pair5, ]

# GLMM FN
GLMM_FN <- pairs_conclusion %>% filter(truth & ancom & zoib & !glmm)
example_pair6 <- c(taxa_pairs$T1[13105],
                   taxa_pairs$T2[13105])
example_counts6 <- counts[example_pair6, ] %>% t() %>% as.data.frame()
example_counts6$Group <- rep(c(0, 1), each=50)
example_counts6$Group <- as.factor(example_counts5$Group)

plot_GLMM_FN <- biplot_viz(example_counts6, example_pair6)
plot_GLMM_FN_combined <- plot_grid(plot_GLMM_FN$countplot,
          plot_GLMM_FN$propplot,
          plot_GLMM_FN$logratioplot,
          ncol=2)

abn_info[example_pair6, ]

ancom_feedback <- ANCOM_result[GLMM_FN$ID, ]
ZOIB_feedback <- ZOIB_result[GLMM_FN$ID, ]
GLMM_feedback <- GLMM_result[GLMM_FN$ID, ]



