# compare GLMdisp and ANCOM in terms of


rm(list=ls())
data_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data'
GLM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/glmdisp_result'
ANCOM_folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/ancom_result'

replicate = 1

performance_summary1 <- function(preddf){
  performance <- table(preddf$truth, preddf$Decision)
  TN <- sum(!preddf$truth & !preddf$Decision)
  FP <- sum(!preddf$truth & preddf$Decision)
  FN <- sum(preddf$truth & !preddf$Decision)
  TP <- sum(preddf$truth & preddf$Decision)
  FDR <- FP/(FP+TP)
  power <- TP/(TP+FN)
  return(list(power=power, FDR=FDR))
}

performance_summary2 <- function(preddf){
  performance <- table(preddf$truth, preddf$prediction)
  TN <- sum(!preddf$truth & !preddf$prediction)
  FP <- sum(!preddf$truth & preddf$prediction)
  FN <- sum(preddf$truth & !preddf$prediction)
  TP <- sum(preddf$truth & preddf$prediction)
  FDR <- FP/(FP+TP)
  power <- TP/(TP+FN)
  return(list(power=power, FDR=FDR))
}

perf_df <- data.frame(Replicate = seq(1, 100),
                      GLM_cutoff_FDR = 0, GLM_cutoff_power=0,
                      ANCOM_cutoff_FDR=0, ANCOM_cutoff_power=0,
                      GLM_SBM_FDR = 0, GLM_SBM_power = 0,
                      ANCOM_SBM_FDR = 0, ANCOM_SBM_power = 0)

for (j in 1:100){
  GLM_cutoff_prediction <- read.csv(file.path(GLM_folder, 'prop_cutoff',
                                              sprintf('glmdisp_decision_%d.csv', j)))
  GLM_cutoff_prediction$truth <- GLM_cutoff_prediction$effect.size != 1
  GLM_cutoff_performance <- performance_summary1(GLM_cutoff_prediction)
  perf_df$GLM_cutoff_FDR[j] <- GLM_cutoff_performance$FDR
  perf_df$GLM_cutoff_power[j] <- GLM_cutoff_performance$power

  GLM_SBM_prediction <- read.csv(file.path(GLM_folder, 'SBM',
                                              sprintf('glmdisp_SBM_decision_%d.csv', j)))
  GLM_SBM_prediction$truth <- GLM_SBM_prediction$effect.size != 1
  GLM_SBM_performance <- performance_summary2(GLM_SBM_prediction)
  perf_df$GLM_SBM_FDR[j] <- GLM_SBM_performance$FDR
  perf_df$GLM_SBM_power[j] <- GLM_SBM_performance$power

  ANCOM_cutoff_prediction <- read.csv(file.path(ANCOM_folder, 'prop_cutoff',
                                                sprintf('ancom_decision_%d.csv', j)))
  ANCOM_cutoff_prediction$truth <- ANCOM_cutoff_prediction$effect.size != 1
  ANCOM_cutoff_performance <- performance_summary1(ANCOM_cutoff_prediction)
  perf_df$ANCOM_cutoff_FDR[j] <- ANCOM_cutoff_performance$FDR
  perf_df$ANCOM_cutoff_power[j] <- ANCOM_cutoff_performance$power

  ANCOM_SBM_prediction <- read.csv(file.path(ANCOM_folder, 'SBM',
                                                sprintf('ancom_SBM_decision_%d.csv', j)))
  ANCOM_SBM_prediction$truth <- ANCOM_SBM_prediction$effect.size != 1
  ANCOM_SBM_performance <- performance_summary2(ANCOM_SBM_prediction)
  perf_df$ANCOM_SBM_FDR[j] <- ANCOM_SBM_performance$FDR
  perf_df$ANCOM_SBM_power[j] <- ANCOM_SBM_performance$power
}

sum(perf_df$GLM_FDR < perf_df$ANCOM_FDR)
sum(perf_df$GLM_power > perf_df$ANCOM_power)




FDR_df <- data.frame(FDR = c(perf_df$GLM_FDR, perf_df$ANCOM_FDR),
                     Method = rep(c("GLMDisp", "ANCOM"), each=100))
FDR_comparison_plot <- ggplot(FDR_df, aes(x=Method, y=FDR)) + geom_violin() +
  geom_boxplot(width=0.02) + ylim(0, 0.60) +
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  xlab("DiffRatio Method") + ylab("False Discovery Rate") +
  theme_bw() + coord_flip()




power_df <- data.frame(power = c(perf_df$GLM_power, perf_df$ANCOM_power),
                       Method = rep(c("GLMDisp", "ANCOM"), each=100))
power_comparison_plot <- ggplot(power_df, aes(x=Method, y=power)) + geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("DiffRatio Method") + ylab("Power") +
  theme_bw() + coord_flip()

library(cowplot)

performance_plot <- plot_grid(FDR_comparison_plot, power_comparison_plot,
                              nrow=1)

# check an example

example_data <- readRDS(file.path(data_folder,
                                  sprintf("simulated_data_%d.rds", 9)))
example_abn <- example_data$mean.eco.abn
example_abn$Taxa <- row.names(example_abn)


example_prediction <- read.csv(file.path(GLM_folder, sprintf('glmdisp_decision_%d.csv', 9)))

positive_decisions <- example_prediction %>% filter(Decision)
min(positive_decisions$SigProp)

FPs <- example_prediction %>% filter(effect.size==1 & Decision) %>%
  inner_join(example_abn, by="Taxa")

TPs <- example_prediction %>% filter(effect.size!=1 & Decision) %>%
  inner_join(example_abn, by="Taxa")

FNs <- example_prediction %>% filter(effect.size!=1 & !Decision) %>%
  inner_join(example_abn, by="Taxa")

TNs <- example_prediction %>% filter(effect.size==1 & !Decision) %>%
  inner_join(example_abn, by="Taxa")
