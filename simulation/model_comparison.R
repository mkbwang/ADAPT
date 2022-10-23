
rm(list=ls())
directory <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'

ancom_performance <- read.csv(file.path(directory, 'ancom_result', 'ancom_performance.csv'))
glmm_performance <- read.csv(file.path(directory, 'glmm_result', 'glmm_performance.csv'))
ZOIB_performance <- read.csv(file.path(directory, 'ZOIB_result', 'ZOIB_performance.csv'))


perform_summary <- function(performance_df){
  pairFDR <- performance_df$pairFP / (performance_df$pairFP + performance_df$pairTP)
  pairpower <- performance_df$pairTP/(performance_df$pairFN + performance_df$pairTP)
  taxapower <- performance_df$taxaTP/(performance_df$taxaTP+performance_df$taxaFN)
  summary_df <- data.frame(pairFDR=pairFDR,
                           pairpower=pairpower,
                           taxapower=taxapower)
  return(summary_df)
}

ancom_summary <- perform_summary(ancom_performance)
ancom_summary$method <- "Wilcoxon"
glmm_summary <- perform_summary(glmm_performance)
glmm_summary$method <- "GLMM"
ZOIB_summary <- perform_summary(ZOIB_performance)
ZOIB_summary$method <- "ZOIB"


combined_summary <- rbind(ancom_summary, glmm_summary, ZOIB_summary)


library(ggplot2)

pairFDR_plot <- ggplot(combined_summary, aes(x=method, y=pairFDR)) +
  geom_boxplot(width=0.15)+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  xlab("Pairwise Testing Method") + ylab("False Discovery Rate of Differential Pairs") +
  coord_flip() + theme_bw()


pairpower_plot <- ggplot(combined_summary, aes(x=method, y=pairpower)) +
  geom_boxplot(width=0.15)+
  xlab("Pairwise Testing Method") + ylab("Power of Detecting Differential Pairs") +
  coord_flip() + theme_bw()


library(cowplot)

pair_performance_plots <- plot_grid(pairFDR_plot, pairpower_plot,
                                    align="v", ncol=1)


# how to cut off to identify differential taxa?
taxapower_plot <- ggplot(combined_summary, aes(x=method, y=taxapower))+
  geom_boxplot(width=0.15)+
  xlab("Pairwise Testing Method") + ylab("Power of Taxa Detection") +
  coord_flip() + theme_bw()
