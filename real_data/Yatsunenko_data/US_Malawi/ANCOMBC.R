library(ANCOMBC)
library(phyloseq)


rm(list=ls())

US_MA_infant_phyobj <- readRDS("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_phyloseq.rds")
US_MA_infant_phyobj <- US_MA_infant_phyobj[-c(4, 10), ]

sample_info <- sample_data(US_MA_infant_phyobj)
otu_mat <- otu_table(US_MA_infant_phyobj)
otu_mat <- otu_mat[-c(4, 10), ]

US_MA_infant_phyobj_subset <- phyloseq(otu_mat, sample_info)

# ancombc
ancombc_result <- ancombc(US_MA_infant_phyobj_subset, formula = 'country', p_adj_method='BH',
                          group='country', struc_zero=FALSE)

decision <- ancombc_result$res$diff_abn


ancombc_summary <- data.frame(Taxa_Name = row.names(decision),
                              Pval = ancombc_result$res$q_val$countryUS,
                              DA = ancombc_result$res$diff_abn$countryUS)

write.csv(ancombc_summary, "real_data/Yatsunenko_data/US_Malawi/US_MA_infant_ANCOMBC_result.csv",
          row.names = FALSE)


