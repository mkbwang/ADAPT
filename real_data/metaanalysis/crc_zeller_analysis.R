rm(list=ls())
library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(DescTools)
source('POLDA/POLDA.R')

zeller_data <- readRDS('real_data/metaanalysis/crc_zeller_results/crc_zeller_phyobj.rds')


# filter out taxa whose prevalence is smaller than 0.05
zeller_filtered_data <- filter_taxa(zeller_data, 
                                    function(x) mean(x != 0) > 0.10 , TRUE)
zeller_metadata <- sample_data(zeller_filtered_data)
zeller_count_mat <- otu_table(zeller_filtered_data)@.Data
zeller_taxonomy_table <- tax_table(zeller_filtered_data)
taxa_with_g <- zeller_taxonomy_table[, 6] !="g__"


# filter out OTU with no names on the genus level
zeller_filtered_data <- subset_taxa(zeller_filtered_data, taxa_with_g)
prevalence <- rowMeans(otu_table(zeller_filtered_data) != 0)
zeller_taxonomy_table <- tax_table(zeller_filtered_data) |> as.matrix()
zeller_metadata <- as(sample_data(zeller_filtered_data), "data.frame")
zeller_metadata$crc <- zeller_metadata$DiseaseState != "H"
sample_data(zeller_filtered_data) <- zeller_metadata
zeller_count_mat <- otu_table(zeller_filtered_data)@.Data

relabd_result <- relabd_GLM(zeller_count_mat, zeller_metadata, "crc")
median_effect <- median(relabd_result$effect)
effect_rank <- rank(abs(relabd_result$effect - median(relabd_result$effect)))
pvals <- relabd_result$pval
AndersonDarlingTest(pvals, null="punif")

subset_taxa <- which(effect_rank < nrow(relabd_result)/2)
taxa_names <- rownames(zeller_count_mat)
zeller_count_mat_subset <- zeller_count_mat[subset_taxa, ]
relabd_result_subset <- relabd_GLM(zeller_count_mat_subset, zeller_metadata, "crc")
subset_pvals <- relabd_result_subset$pval
AndersonDarlingTest(subset_pvals, null="punif")


other_taxa_result <- reference_GLM(zeller_count_mat, zeller_metadata, "crc",
                                   reftaxa=taxa_names[subset_taxa])

final_result <- rbind(relabd_result_subset, other_taxa_result)

# q value adjustment
qobj <- qvalue(final_result$pval)

# BH adjustment
BH_pval <- p.adjust(other_taxa_result$pval, method="BH")

# BUM estimation
library(BioNet)
BUM_fit <- fitBumModel(final_result$pval)
lambda_est <- BUM_fit$lambda
a_est <- BUM_fit$a 
pi_ub <- lambda_est + (1-lambda_est) * a_est

cutoff_calculation <- function(a_est, lambda_est, alpha){
  pi_ub <- lambda_est + (1-lambda_est) * a_est
  tau <- ((pi_ub - alpha*lambda_est)/(alpha * (1-lambda_est)))^(1/(a_est - 1))
  return(tau)
}

# 
# POLDA analysis
polda_result <- polda(otu_table = zeller_count_mat,
                       metadata = zeller_metadata,
                       covar="crc")
polda_df <- polda_result$P_Value |> as.data.frame()
polda_DA_taxa <- polda_result$DA_taxa


# plot taxa prevalence
prevalence <- rowMeans(zeller_count_mat != 0)
prevalence_df <- data.frame(Prevalence=prevalence)
hist_prevalence <- ggplot(prevalence_df, aes(x=Prevalence)) +
  geom_histogram(color="black", fill="white", bins=20) +
  xlab('') + ggtitle("Genus Prevalence") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


# plot sequencing depths
library_size <- colSums(zeller_count_mat)
library_size_df <- data.frame(Library_Size=library_size)
hist_library_size <- ggplot(library_size_df, aes(x=Library_Size))+
  geom_histogram(color="black", fill="white", bins=20) +
  scale_x_log10(breaks = c(30000, 50000, 1e5, 2e5, 3e5)) +
  xlab('')+
  ggtitle("Sequencing Depths") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
basic_stats_plot <- plot_grid(hist_prevalence, hist_library_size, nrow=1)

# separated steps
zeller_glm_alltaxa <- relabd_GLM(zeller_count_mat,
                               zeller_metadata,
                               covar="crc")


# plot estimated tau for the full model
tau_hats <- zeller_glm_alltaxa$effect
median_tau_hat <- median(tau_hats)
Taxa_info_df <- data.frame(Taxa = zeller_glm_alltaxa$Taxa,
                           Tau = tau_hats)
Taxa_info_df$Prevalence <- prevalence[Taxa_info_df$Taxa]
hist_tau <- ggplot(Taxa_info_df, aes(x=Tau)) +
  geom_histogram(color="black", fill="white", bins=20) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(xintercept=median_tau_hat, linetype="dashed", color="red")+
  xlab('') + ggtitle("Estimated Tau") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


# Fit ANCOMBC
ancombc_result <- ancombc(zeller_filtered_data, formula = 'crc', p_adj_method='BH',
                          group='crc', struc_zero=TRUE,
                          zero_cut=0.95)
ancombc_call <- ancombc_result$res |> as.data.frame()
colnames(ancombc_call) <- c("beta", "se", "W", "p_val", "q_val", "diff_abn")
ancombc_DA_taxa <- rownames(ancombc_call)[ancombc_call$diff_abn]



# Fit Aldex2
aldex_result <- aldex(zeller_count_mat, zeller_metadata$crc, test='t')
aldex_DA_taxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]

zeller_taxonomy_df <- as(zeller_taxonomy_table, 'matrix') |> as.data.frame()
zeller_taxonomy_df$polda_effect <- NA
zeller_taxonomy_df[polda_df$Taxon, 'polda_effect'] <- polda_df$effect
zeller_taxonomy_df$polda_rawp <- NA
zeller_taxonomy_df[polda_df$Taxon, 'polda_rawp'] <- polda_df$pval
zeller_taxonomy_df$polda_adjustedp <- NA
zeller_taxonomy_df[polda_df$Taxon, 'polda_adjustedp'] <- polda_df$pval_adjust
zeller_taxonomy_df$ancombc_rawp <- ancombc_call$p_val
zeller_taxonomy_df$ancombc_adjustedp <- ancombc_call$q_val
zeller_taxonomy_df$aldex_rawp <- aldex_result$wi.ep
zeller_taxonomy_df$aldex_adjustedp <- aldex_result$wi.eBH



write.csv(zeller_taxonomy_df, 'real_data/metaanalysis/zeller_performance.csv')


