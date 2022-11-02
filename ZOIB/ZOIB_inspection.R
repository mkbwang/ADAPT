library(zoib)
library(dplyr)
#--------------pvalue inspection on local machine---------------#

data_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/RMD/CAARS_data'
model_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/RMD/CAARS_Model_Summary'

counts_table <- read.csv(file.path(data_folder, 'filtered_count_order_level.csv'),
                         row.names = 1)
metadata <- read.csv(file.path(data_folder, 'filtered_metadata_order_level.csv'))

pvals <- read.csv(file.path(model_folder, 'Pvals_order_level.csv'))
pvals_zoib <- read.csv(file.path(model_folder, 'Pvals_order_level_zoib.csv'))

zoib_estimation <- function(j){
  t1 <- pvals$Taxa1[j]
  t2 <- pvals$Taxa2[j]

  counts_t1 <- counts_table[, t1]
  counts_t2 <- counts_table[, t2]

  t1pt2n <- (counts_t1 > 0) & (counts_t2 == 0)
  t1nt2p <- (counts_t1 == 0) & (counts_t2 > 0)
  t1prop <- counts_t1 / (counts_t1 + counts_t2)

  df <- cbind(metadata$asthma, t1pt2n, t1nt2p, t1prop, counts_t1, counts_t2) %>% as.data.frame()
  colnames(df) <- c("Asthma", "T1PT2N", 'T1NT2P', 'T1prop', 'count_t1', 'count_t2')

  # first check if there is difference in zero inflation
  zinfmodel <- glm(T1NT2P ~ Asthma, data=df,
                   family = binomial(link = "logit"))
  zinfpval <- summary(zinfmodel)$coefficients[2,4]
  # pvals_table$infl0[j] <- zinfpval

  df_filtered <- df %>% filter(!T1NT2P)
  oinfmodel <- glm(T1PT2N ~ Asthma, data=df_filtered,
                   family = binomial(link = "logit"))
  oinfpval <- summary(oinfmodel)$coefficients[2,4]
  # pvals_table$infl1[j] <- oinfpval

  df_dp <- df %>% filter(count_t1 > 0 & count_t2 > 0)
  betapval <- NA
  # if (nrow(df_dp) > 1){
  #   invisible(betamodel <- zoib(T1prop ~ Asthma|Asthma,
  #                               data=df_dp, random=0, zero.inflation = FALSE,
  #                               one.inflation = FALSE, joint=FALSE,
  #                               n.iter=8000, n.thin=15, n.burn=2000))
  #   betareg_summary <- summary(betamodel$coeff)$statistics
  #   beta_test_statistic <- -abs(betareg_summary[2,1] / betareg_summary[2,2])
  #   betapval <- pnorm(beta_test_statistic)
  #   # pvals_table$betareg[j] <- betapval
  # }
  result = list(Taxa1=t1, Taxa2=t2, Inf0 = zinfpval, Inf1=oinfpval,
                betapval=betapval)
  return(as.data.frame(result))
}




pvals$zinf <- NA
pvals$oinf <- NA
pvals$betareg <- NA
pvals$zinf[pvals_zoib$Rownum] <- pvals_zoib$Inf0
pvals$oinf[pvals_zoib$Rownum] <- pvals_zoib$Inf1
pvals$betareg[pvals_zoib$Rownum] <- pvals_zoib$betapval

## for some reason the results generated from cluster shifted by one row

result_last <- zoib_estimation(378)
# result_second_last <- zoib_estimation(377)
# result_third_last <- zoib_estimation(376)

pvals$zinf[1:377] <- pvals$zinf[2:378]
pvals$oinf[1:377] <- pvals$oinf[2:378]
pvals$betareg[1:377] <- pvals$betareg[2:378]
pvals$zinf[378] <- result_last$Inf0
pvals$oinf[378] <- result_last$Inf1
pvals$betareg[378] <- result_last$betapval

# check the pairs that failed, for some reason only one of them was actually broken
failed_pair <- which(is.na(pvals$zinf))

for (rowid in failed_pair){
  result <- zoib_estimation(rowid)
  pvals$zinf[rowid] <- result$Inf0
  pvals$oinf[rowid] <- result$Inf1
  pvals$betareg[rowid] <- result$betapval
}


# compare the adjusted p values between ANCOM, GLMM and ZOIB

ANCOM_pvals <- pvals$wilx
ANCOM_pvals_adjusted <- p.adjust(ANCOM_pvals, method="BH")
GLMM_pvals <- pvals$GLMM_GQ_BOBYQA
GLMM_pvals_adjusted <- p.adjust(GLMM_pvals, method="BH")

zinf_pval <- pvals$zinf
zinf_pval_subset <- zinf_pval[zinf_pval < 1]
zinf_pval_adjusted <- p.adjust(zinf_pval_subset, method="BH")
oinf_pval <- pvals$oinf
oinf_pval_subset <- oinf_pval[oinf_pval < 1]
oinf_pval_adjusted <- p.adjust(oinf_pval_subset, method="BH")
beta_pval <- pvals$betareg
beta_pval_adjusted <- p.adjust(beta_pval, method="BH")

library(ggplot2)
tpairs_significant <- which(beta_pval_adjusted < 0.05)
plot_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/RMD/CAARS_plots'

biplot_viz <- function(tpair){
  subset_data <- cbind(metadata$asthma, counts_table[, tpair]) %>%
    as.data.frame()
  colnames(subset_data) <- c("Asthma", tpair)
  subset_data$ratios = subset_data[, 2] / (subset_data[, 2] + subset_data[, 3])
  subset_data$Asthma <- as.factor(subset_data$Asthma)

  plot1_name <- sprintf('%s_%s_biplot.pdf', tpair[1], tpair[2])
  count_plot <- ggplot(subset_data, aes_string(x=tpair[1], y=tpair[2], color='Asthma')) +
    geom_jitter(alpha=0.6) +
    theme_bw()
  ggsave(file.path(plot_folder, 'biplots_significant', plot1_name),
         plot=count_plot,
         width=12, height=8, units="cm")

  plot2_name <- sprintf('%s_%s_Ratios.pdf', tpair[1], tpair[2])
  ratio_plot <- ggplot(subset_data, aes(x=ratios, color=Asthma)) +
    geom_histogram(fill="white", binwidth = 0.05)+
    xlab(sprintf('%s/(%s+%s)', tpair[1], tpair[1], tpair[2])) + theme_bw()
  ggsave(file.path(plot_folder, 'ratios_significant', plot2_name),
         plot=ratio_plot,
         width=12, height=8, units="cm")
}

for (i in tpairs_significant){
  selected_t1 <- pvals$Taxa1[i]
  selected_t2 <- pvals$Taxa2[i]
  biplot_viz(c(selected_t1, selected_t2))
}

View(pvals[beta_pval_adjusted < 0.05, ])
write.csv(pvals, file.path(model_folder, 'Pvals_order_level.csv'), row.names = FALSE)

