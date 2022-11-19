
rm(list=ls())
directory <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'

ancom_performance <- read.csv(file.path(directory, 'ancom_result', 'ancom_performance.csv'))
glmdisp_performance <- read.csv(file.path(directory, 'glmdisp_result', 'glmdisp_performance.csv'))
ZOIB_performance <- read.csv(file.path(directory, 'ZOIB_result', 'ZOIB_performance.csv'))


perform_summary <- function(performance_df){
  pairFDR <- performance_df$pairFP / (performance_df$pairFP + performance_df$pairTP)
  pairpower <- performance_df$pairTP/(performance_df$pairFN + performance_df$pairTP)
  taxaFDR <- performance_df$taxaFP / (performance_df$taxaFP + performance_df$taxaTP)
  taxapower <- performance_df$taxaTP/(performance_df$taxaTP+performance_df$taxaFN)
  summary_df <- data.frame(pairFDR=pairFDR,
                           pairpower=pairpower,
                           taxaFDR=taxaFDR,
                           taxapower=taxapower)
  return(summary_df)
}

ancom_summary <- perform_summary(ancom_performance)
ancom_summary$method <- "Wilcoxon"
glmdisp_summary <- perform_summary(glmdisp_performance)
glmdisp_summary$method <- "GLMdisp"
ZOIB_summary <- perform_summary(ZOIB_performance)
ZOIB_summary$method <- "ZOIB"

deseq2_summary <- read.csv(file.path(directory, 'DESEQ2_performance.csv')) %>%
  select(FDR, Power)
colnames(deseq2_summary) <- c("taxaFDR", "taxapower")
deseq2_summary$method <- "DESeq2"


combined_summary <- rbind(ancom_summary, glmdisp_summary, ZOIB_summary)


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


taxa_performance <- rbind(combined_summary %>% select(taxaFDR, taxapower, method),
                          deseq2_summary)

# how to cut off to identify differential taxa?
taxaFDR_plot <- ggplot(taxa_performance, aes(x=method, y=taxaFDR))+
  geom_boxplot(width=0.15)+geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  xlab("Pairwise Testing Method") + ylab("FDR of Differential Taxa Detection") +
  coord_flip() + theme_bw()


taxapower_plot <- ggplot(taxa_performance, aes(x=method, y=taxapower))+
  geom_boxplot(width=0.15)+
  xlab("Pairwise Testing Method") + ylab("Power of Differential Taxa Detection") +
  coord_flip() + theme_bw()


taxa_performance_plots <- plot_grid(taxaFDR_plot, taxapower_plot,
                                    align="v", ncol=1)
