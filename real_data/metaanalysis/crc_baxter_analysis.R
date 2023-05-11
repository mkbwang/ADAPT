rm(list=ls())
library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
source('POLDA/POLDA.R')
baxter_data <- readRDS('real_data/metaanalysis/crc_baxter_results/crc_baxter_phyobj.rds')


# filter out taxa whose prevalence is smaller than 0.05
baxter_filtered_data <- filter_taxa(baxter_data, 
                                      function(x) mean(x != 0) > 0.05 , TRUE)
baxter_metadata <- sample_data(baxter_filtered_data)
prevalence <- rowMeans(otu_table(baxter_filtered_data) != 0)
baxter_count_mat <- otu_table(baxter_filtered_data)@.Data
baxter_taxonomy_table <- tax_table(baxter_filtered_data)
taxa_with_g <- baxter_taxonomy_table[, 6] !="g__"

# filter out OTU with no names on the genus level
baxter_filtered_data <- subset_taxa(baxter_filtered_data, taxa_with_g)
prevalence <- rowMeans(otu_table(baxter_filtered_data) != 0)
baxter_taxonomy_table <- tax_table(baxter_filtered_data) |> as.matrix()
baxter_metadata <- as(sample_data(baxter_filtered_data), "data.frame")
baxter_metadata$crc <- baxter_metadata$DiseaseState != "H"
sample_data(baxter_filtered_data) <- baxter_metadata
baxter_count_mat <- otu_table(baxter_filtered_data)@.Data
zero_counts <- rowSums(baxter_count_mat == 0)

# prevalence histogram
prevalence <- rowMeans(baxter_count_mat != 0)
prevalence_df <- data.frame(Prevalence=prevalence)
hist_prevalence <- ggplot(prevalence_df, aes(x=Prevalence)) +
  geom_histogram(color="black", fill="white", bins=20) +
  xlab('') + ggtitle("Genus Prevalence") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
# sequence depths histogram
library_size <- colSums(baxter_count_mat)
library_size_df <- data.frame(Library_Size=library_size)
hist_library_size <- ggplot(library_size_df, aes(x=Library_Size))+
  geom_histogram(color="black", fill="white", bins=20) +
  scale_x_log10(breaks = c(1000, 5000, 1e4, 5e4, 1e5)) +
  xlab('')+
  ggtitle("Sequencing Depths") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
basic_stats_plot <- plot_grid(hist_prevalence, hist_library_size, nrow=1)


# three steps for polda
baxter_pairglm <- pairwise_GLM(baxter_count_mat,
                               baxter_metadata,
                               covar="crc")

reftaxa_search <- backward_selection(baxter_count_mat,
                                     baxter_pairglm,
                                     start="median")

# plot estimated tau hats for all the taxa
tau_hats <- reftaxa_search$full_tau_hat
median_tau_hat <- median(tau_hats)
Taxa_info_df <- data.frame(Taxa = names(tau_hats),
                           Tau = tau_hats)

Taxa_info_df$Prevalence <- prevalence[Taxa_info_df$Taxa]
hist_tau <- ggplot(Taxa_info_df, aes(x=Tau)) +
  geom_histogram(color="black", fill="white", bins=25) +
  scale_x_continuous(breaks=seq(-3, 3, 0.5))+
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(xintercept=median_tau_hat, linetype="dashed", color="red")+
  xlab('') + ggtitle("Estimated Tau") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


# plot prevalence over the tau hat 
Taxa_subset <- Taxa_info_df %>% filter(abs(Tau - median_tau_hat) < 0.08)
Taxa_subset <- Taxa_subset %>% arrange(Tau)
prev_tau_plot <- ggplot(Taxa_subset, aes(x=Tau, y=Prevalence)) +
  geom_point(size=1.5) + 
  geom_vline(xintercept=median_tau_hat, linetype="dashed", color="red")+
  scale_x_continuous(breaks=seq(-0.17, 0.01, 0.02))+
  geom_label_repel(aes(label=Taxa), size=2.5) + xlab("Tau") + ylab("Prevalence")+
  theme_bw()


polda_result <- reference_GLM(baxter_count_mat, baxter_metadata, 'crc',
                              reftaxa=reftaxa_search$reftaxa)
# Fit ANCOMBC
ancombc_result <- ancombc(baxter_filtered_data, formula = 'crc', p_adj_method='BH',
                          group='crc', struc_zero=TRUE,
                          zero_cut=0.95)
ancombc_call <- ancombc_result$res |> as.data.frame()
colnames(ancombc_call) <- c("beta", "se", "W", "p_val", "q_val", "diff_abn")
ancombc_DA_taxa <- rownames(ancombc_call)[ancombc_call$crcTRUE]


# Fit Aldex2
aldex_result <- aldex(baxter_count_mat, baxter_metadata$crc, test='t')
aldex_DA_taxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]

baxter_taxonomy_df <- as(baxter_taxonomy_table, 'matrix') |> as.data.frame()
baxter_taxonomy_df$polda_effect <- NA
baxter_taxonomy_df[polda_result$Taxon, 'polda_effect'] <- polda_result$effect
baxter_taxonomy_df$polda_rawp <- NA
baxter_taxonomy_df[polda_result$Taxon, 'polda_rawp'] <- polda_result$pval
baxter_taxonomy_df$polda_adjustedp <- NA
baxter_taxonomy_df[polda_result$Taxon, 'polda_adjustedp'] <- polda_result$pval_adjust
baxter_taxonomy_df$ancombc_rawp <- ancombc_call$p_val
baxter_taxonomy_df$ancombc_adjustedp <- ancombc_call$q_val
baxter_taxonomy_df$aldex_rawp <- aldex_result$wi.ep
baxter_taxonomy_df$aldex_adjustedp <- aldex_result$wi.eBH


write.csv(baxter_taxonomy_df, 'real_data/metaanalysis/baxter_performance.csv')

