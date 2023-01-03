# empirical way of p value cutoff

library(dplyr)
data_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data'

GLM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/glmdisp_result'
ANCOM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/ancom_result'

# decide the proportion cutoff
threshold_selection <- function(significant_proportions){
  significant_proportions<- significant_proportions %>% arrange(SigProp)
  sorted_proportions_count <- significant_proportions %>% group_by(SigProp) %>%
    summarise(count = n())
  sorted_proportions_count$cumcount <- cumsum(sorted_proportions_count$count)
  sorted_proportions_count$cumcount <- sorted_proportions_count$cumcount / max(sorted_proportions_count$cumcount)
  sorted_proportions_count$criteria <- 1 - sorted_proportions_count$SigProp - sorted_proportions_count$cumcount
  selected_index <- min(which(sorted_proportions_count$criteria <= 0))
  threshold <- sorted_proportions_count$SigProp[selected_index]
  return(list(sorted_proportions_count, threshold))
}

# generate plot
library(ggplot2)
library(cowplot)

gen_plot <- function(cutoff_result, sigprop_df){

  sorted_proportions_count <- cutoff_result[[1]]
  threshold <- cutoff_result[[2]]

  sorted_proportions_plot <- data.frame(Proportion = c(0, rep(sorted_proportions_count$SigProp, each=2), 1),
                                   CumProb = rep(c(0, sorted_proportions_count$cumcount), each=2))
  histogram_cutoff_plot <- ggplot(sigprop_df, aes(x=SigProp)) +
    geom_histogram(color="black", fill="white", binwidth=0.02) +
    scale_x_continuous(breaks=seq(0, 1, 0.1))+
    geom_vline(xintercept=threshold, linetype="dashed", color="red") +
    xlab("Proportion of Differential Pairs") + ylab("Frequency")+
    theme_bw()

  sorted_proportions_plot$Type <- "CumProb"
  diagline <- data.frame(Proportion=c(0, 1),
                         CumProb=c(1, 0),
                         Type=c("Assist", "Assist"))
  plot_df <- rbind(sorted_proportions_plot, diagline)


  derivation_plot <- ggplot(plot_df, aes(x=Proportion, y=CumProb)) +
    geom_line(aes(linetype=Type)) +
    scale_linetype_manual(values=c("twodash", "solid"))+
    theme_bw()+theme(legend.position="none")+xlab("Proportion of Differential Pairs")+
    ylab("Cumulative Probability")

  return(list(derivation = derivation_plot, illustration = histogram_cutoff_plot))
}


library(foreach)
library(doParallel)
library(dplyr)
cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

replicate = 1
for (replicate in 1:100){
  # get the truth
  simulated_data = readRDS(file.path(data_folder,
                                     sprintf("simulated_data_%d.rds", replicate)))
  abn_info <- simulated_data$mean.eco.abn
  truth <- data.frame(Taxa = row.names(abn_info),
                      effect.size = abn_info$effect.size)

  # load GLMDISP result
  GLM_result = read.csv(file.path(GLM_folder,
                                  sprintf("glmdisp_result_%d.csv", replicate)))
  GLM_result$adjusted_pval <- p.adjust(GLM_result$pval, method='BH')

  taxa_names <- sprintf('taxon%d', seq(1, 1000))


  # calculate the significant pair proportion for each taxon
  count_significant_prop <- foreach(j=1:1000, .packages = c("dplyr"), .combine=c) %dopar%{
    selected_taxon <- taxa_names[j]
    relevant_pairs <- GLM_result %>%
      filter(T1 == selected_taxon | T2 == selected_taxon)
    return(mean(relevant_pairs$adjusted_pval < 0.05))
  }


  significant_proportions <- data.frame(Taxa = taxa_names,
                                        SigProp = count_significant_prop)


  cutoff <- threshold_selection(significant_proportions)
  threshold <- cutoff[[2]]

  visualizations <- gen_plot(cutoff, significant_proportions)
  combined_plots <- plot_grid(visualizations$derivation, visualizations$illustration,
                              nrow=1)

  significant_proportions$Decision <- significant_proportions$SigProp >= threshold
  result <- truth %>% full_join(significant_proportions, by="Taxa")
  write.csv(result, file.path(GLM_folder, sprintf('glmdisp_decision_%d.csv', replicate)),
            row.names=FALSE)
  ggsave(file.path(GLM_folder, sprintf('glmdisp_decision_%d.pdf', replicate)),
         combined_plots, width=15, height=8, units="cm")
}


for (replicate in 1:100){
  # get the truth
  simulated_data = readRDS(file.path(data_folder,
                                     sprintf("simulated_data_%d.rds", replicate)))
  abn_info <- simulated_data$mean.eco.abn
  truth <- data.frame(Taxa = row.names(abn_info),
                      effect.size = abn_info$effect.size)

  # load GLMDISP result
  ANCOM_result = read.csv(file.path(ANCOM_folder,
                                  sprintf("ancom_result_%d.csv", replicate)))
  ANCOM_result$adjusted_pval <- p.adjust(ANCOM_result$pval, method='BH')

  taxa_names <- sprintf('taxon%d', seq(1, 1000))


  # calculate the significant pair proportion for each taxon
  count_significant_prop <- foreach(j=1:1000, .packages = c("dplyr"), .combine=c) %dopar%{
    selected_taxon <- taxa_names[j]
    relevant_pairs <- ANCOM_result %>%
      filter(T1 == selected_taxon | T2 == selected_taxon)
    return(mean(relevant_pairs$adjusted_pval < 0.05))
  }


  significant_proportions <- data.frame(Taxa = taxa_names,
                                        SigProp = count_significant_prop)


  cutoff <- threshold_selection(significant_proportions)
  threshold <- cutoff[[2]]

  visualizations <- gen_plot(cutoff, significant_proportions)
  combined_plots <- plot_grid(visualizations$derivation, visualizations$illustration,
                              nrow=1)

  significant_proportions$Decision <- significant_proportions$SigProp >= threshold
  result <- truth %>% full_join(significant_proportions, by="Taxa")
  write.csv(result, file.path(ANCOM_folder, sprintf('ancom_decision_%d.csv', replicate)),
            row.names=FALSE)
  ggsave(file.path(ANCOM_folder, sprintf('ancom_decision_%d.pdf', replicate)),
         combined_plots, width=15, height=8, units="cm")
}


# ANCOM_result <- read.csv(file.path(ANCOM_folder,
#                                    sprintf("ancom_result_%d.csv", replicate)))
# ANCOM_result$adjusted_pval <- p.adjust(ANCOM_result$pval, method='BH')



