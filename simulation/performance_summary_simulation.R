rm(list=ls())
# performance_summary
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
data_folder <- '/home/wangmk/MDAWG/POLDA/simulation/data'
aldex_folder <- '/home/wangmk/MDAWG/POLDA/simulation/Aldex2'
ANCOM_folder <- '/home/wangmk/MDAWG/POLDA/simulation/ANCOM'
ANCOMBC_folder <- '/home/wangmk/MDAWG/POLDA/simulation/ANCOMBC'
# DACOMP_folder <- '/home/wangmk/MDAWG/POLDA/simulation/DACOMP'
POLDA_folder <- "/home/wangmk/MDAWG/POLDA/simulation/POLDA"


# polda result summary
polda_performance <- data.frame(ID=seq(1, 100),
                                AUROC=rep(0, 100),
                                FDR=rep(0, 100),
                                Power=rep(0, 100))
for (j in 1:100){
  polda_result <- read.csv(file.path(POLDA_folder, 
                                     sprintf('polda_simulation_result_%d.csv', j)))
  ## for those reference taxa, I put a placeholder of 1 as their p values
  polda_result$pval_adjust[polda_result$RefTaxa] <- 1
  polda_result_subset <- polda_result[!polda_result$RefTaxa, ]
  ## if there are still taxa whose p values are NA, then they are structural zeros
  polda_result$pval_adjust[is.na(polda_result$pval_adjust)] <- 0
  polda_roc <- roc(polda_result_subset$logfold!=0, 1-polda_result_subset$pval_adjust)
  polda_performance$AUROC[j] <- polda_roc$auc
  DA_truth <- polda_result$logfold != 0
  polda_contingency <- table(DA_truth, polda_result$pval_adjust < 0.05)
  if (ncol(polda_contingency) == 2){
    polda_performance$FDR[j] <- polda_contingency[1, 2] / sum(polda_contingency[, 2])
    polda_performance$Power[j] <- polda_contingency[2, 2] / sum(polda_contingency[2, ]) 
  }
}
polda_performance$Method <- 'POLDA'


# ancombc result summary
ancombc_performance <- data.frame(ID=seq(1, 100),
                                  AUROC=rep(0, 100),
                                  FDR=rep(0, 100),
                                  Power=rep(0, 100))

for (j in 1:100){
  ancombc_result <- read.csv(file.path(ANCOMBC_folder, sprintf('ancombc_simulation_result_%d.csv', j)))
  ## If P value is NA, that means structural zeros. I put a zero p value as place holder
  ancombc_result$Pval_Adjust[is.na(ancombc_result$Pval_Adjust)] <- 0
  DA_truth <- ancombc_result$logfold !=0
  ancombc_roc <- roc(DA_truth, 1-ancombc_result$Pval_Adjust)
  ancombc_performance$AUROC[j] <- ancombc_roc$auc
  ancombc_contingency <- table(DA_truth, ancombc_result$Pval_Adjust < 0.05)
  ancombc_performance$FDR[j] <- ancombc_contingency[1,2]/sum(ancombc_contingency[, 2])
  ancombc_performance$Power[j] <- ancombc_contingency[2,2]/sum(ancombc_contingency[2,])
}
ancombc_performance$Method="ANCOM-BC"

# Aldex2 summary
aldex2_performance <- data.frame(ID=seq(1, 100),
                                 AUROC=rep(0, 100),
                                 FDR=rep(0, 100),
                                 Power=rep(0, 100))

for (j in 1:100){
  aldex2_result <- read.csv(file.path(aldex_folder, sprintf('aldex_simulation_result_%d.csv', j)))
  DA_truth <- aldex2_result$logfold != 0
  aldex2_roc <- roc(DA_truth, 1-aldex2_result$DA_Pval)
  aldex2_performance$AUROC[j] <- aldex2_roc$auc
  aldex2_contingency <- table(DA_truth, aldex2_result$DA_Pval < 0.05)
  if (ncol(aldex2_contingency) == 2){
    aldex2_performance$FDR[j] <- aldex2_contingency[1,2]/sum(aldex2_contingency[,2])
    aldex2_performance$Power[j] <- aldex2_contingency[2, 2] / sum(aldex2_contingency[2, ]) 
  }
}
aldex2_performance$Method <- "Aldex2"

#ANCOM Summary
ancom_performance <- data.frame(ID=seq(1, 100),
                                AUROC=rep(0, 100),
                                FDR=rep(0, 100),
                                Power=rep(0, 100))

for (j in 1:100){
  ancom_result <- read.csv(file.path(ANCOM_folder, sprintf('ancom_simulation_result_%d.csv', j)))
  ## if a taxon has W being NA, it means structural zeros
  ancom_result$W[is.na(ancom_result$W)] <- nrow(ancom_result) - 1
  DA_truth <- ancom_result$logfold != 0
  ancom_roc <- roc(DA_truth, ancom_result$W)
  ancom_performance$AUROC[j] <- ancom_roc$auc
}
ancom_performance$Method <- "ANCOM"

auroc_combined <- cbind(aldex2_performance$AUROC,
                        ancom_performance$AUROC,
                        ancombc_performance$AUROC,
                        polda_performance$AUROC) |> as.data.frame()
colnames(auroc_combined) <- c("Aldex2", 'ANCOM', "ANCOMBC", "POLDA")
winner <- apply(auroc_combined, 1, which.max)

auroc_diff_df <- cbind(polda_performance$AUROC - aldex2_performance$AUROC,
                       polda_performance$AUROC - ancombc_performance$AUROC,
                       polda_performance$AUROC - ancom_performance$AUROC) |> as.data.frame()

colnames(auroc_diff_df) <- c("POLDA-ALDEx2", "POLDA-ANCOMBC", "POLDA-ANCOM")
auroc_diff_df_long <- auroc_diff_df %>% pivot_longer(cols=c("POLDA-ALDEx2", "POLDA-ANCOMBC", "POLDA-ANCOM"),
                                          names_to='Pair',
                                          values_to='Difference')



combined_performance <- rbind(aldex2_performance, ancom_performance, ancombc_performance, polda_performance)

library(cowplot)

AUROC_comparison_plot1 <- ggplot(combined_performance, aes(x=Method, y=AUROC)) + geom_violin()+
  geom_boxplot(width=0.02) + xlab("") + ylab("AUROC")+theme_bw()+coord_flip()

AUROC_comparison_plot2 <- ggplot(auroc_diff_df_long, aes(x=Pair, y=Difference)) + geom_violin()+
  geom_boxplot(width=0.02) +geom_hline(yintercept=0, linetype="dashed", color="red") +
  xlab("") + ylab("AUROC Difference")+theme_bw()+coord_flip()

plot_grid(AUROC_comparison_plot1, AUROC_comparison_plot2, ncol=1)

combined_performance_subset <- combined_performance %>% filter(Method != "ANCOM")

FDR_comparison_plot <- ggplot(combined_performance_subset, aes(x=Method, y=FDR)) + geom_violin() +
  geom_boxplot(width=0.02) + ylim(0, 0.4) +
  geom_hline(yintercept=0.05, linetype="dashed", color="red")+
  xlab("") + ylab("False Discovery Rate") +
  theme_bw() + coord_flip()

power_comparison_plot <- ggplot(combined_performance_subset, aes(x=Method, y=Power)) + geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("") + ylab("Power") +
  theme_bw() + coord_flip()
library(cowplot)
plot_grid(FDR_comparison_plot, power_comparison_plot, ncol=1)

# POLDA_taxa_count$Falseref <- POLDA_taxa_count$Falseref > 0
# ggplot(POLDA_taxa_count, aes(x=RefTaxa, y=DiffTaxa, color=Falseref)) +
#   geom_jitter(size=1.5, alpha=0.7, width=0.1) + theme_bw() +
#   scale_color_manual(values=c("gray1", "brown2")) +
#   xlab("Number of Reference Taxa") + ylab("Number of DA Taxa")+
#   theme(legend.position = "None")
