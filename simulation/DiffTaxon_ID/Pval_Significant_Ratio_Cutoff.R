# empirical way of p value cutoff

library(dplyr)
data_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data'

pairGLM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/glmdisp_result'
ANCOM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/ancom_result'

replicate = 1

# get the truth
simulated_data = readRDS(file.path(data_folder,
                                   sprintf("simulated_data_%d.rds", replicate)))
abn_info <- simulated_data$mean.eco.abn
truth <- data.frame(Taxa = row.names(abn_info),
                    effect.size = abn_info$effect.size)

# load GLMDISP result
GLM_result = read.csv(file.path(pairGLM_folder,
                                sprintf("glmdisp_result_%d.csv", replicate)))
GLM_result$adjusted_pval <- p.adjust(GLM_result$pval, method='BH')

taxa_names <- sprintf('taxon%d', seq(1, 1000))
library(foreach)
library(doParallel)
library(dplyr)
cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


# calculate the significant pair proportion for each taxon
count_significant_ratio <- foreach(j=1:1000, .packages = c("dplyr"), .combine=c) %dopar%{
  selected_taxon <- taxa_names[j]
  relevant_pairs <- GLM_result %>%
    filter(T1 == selected_taxon | T2 == selected_taxon)
  return(mean(relevant_pairs$adjusted_pval < 0.05))
}


significant_proportions <- data.frame(Taxa = taxa_names,
                                 SigProp = count_significant_ratio)


# decide the proportion cutoff
threshold_selection <- function(significant_proportions){
  significant_ratios <- significant_ratios %>% arrange(SigRatio)
  sorted_ratios_count <- significant_ratios %>% group_by(SigRatio) %>%
    summarise(count = n())
  sorted_ratios_count$cumcount <- cumsum(sorted_ratios_count$count)
  sorted_ratios_count$cumcount <- sorted_ratios_count$cumcount / max(sorted_ratios_count$cumcount)
  sorted_ratios_count$criteria <- 1 - sorted_ratios_count$SigRatio - sorted_ratios_count$cumcount
  selected_index <- min(which(sorted_ratios_count$criteria <= 0))
  threshold <- sorted_ratios_count$SigRatio[selected_index]
  return(list(sorted_ratios_count, threshold))
}



cutoff <- threshold_selection(significant_proportions)
threshold <- cutoff[[2]]

# generate plot
library(ggplot2)

gen_plot <- function(sorted_ratios_count){

  sorted_ratios_plot <- data.frame(Proportion = c(0, rep(sorted_ratios_count$SigRatio, each=2), 1),
                                   CumProb = rep(c(0, sorted_ratios_count$cumcount), each=2))
  sorted_ratios_plot$Type <- "CumProb"
  diagline <- data.frame(Proportion=c(0, 1),
                         CumProb=c(1, 0),
                         Type=c("Assist", "Assist"))
  plot_df <- rbind(sorted_ratios_plot, diagline)


  illustration <- ggplot(plot_df, aes(x=Proportion, y=CumProb)) +
    geom_line(aes(linetype=Type)) +
    scale_linetype_manual(values=c("twodash", "solid"))+
    theme_bw()+theme(legend.position="none")+xlab("Significant Pair Proportion")+
    ylab("Cumulative Probability")

  return(illustration)
}


illustration_cutoff <- gen_plot(cutoff[[1]])
significant_proportions$Decision <- significant_proportions$SigProp >= threshold
result <- truth %>% full_join(significant_proportions, by="Taxa")
write.csv(result, file.path(pairGLM_folder, sprintf('glmdisp_decision_%d.csv', replicate)),
          row.names=FALSE)
ggsave(file.path(pairGLM_folder, sprintf('glmdisp_decision_%d.pdf', replicate)),
       illustration_cutoff, width=10, height=8, units="cm")


# ANCOM_result <- read.csv(file.path(ANCOM_folder,
#                                    sprintf("ancom_result_%d.csv", replicate)))
# ANCOM_result$adjusted_pval <- p.adjust(ANCOM_result$pval, method='BH')



