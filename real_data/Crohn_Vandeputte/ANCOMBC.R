
library(ANCOMBC)
library(phyloseq)

rm(list=ls())

crohn_phyobj <- readRDS("real_data/Crohn_Vandeputte/Crohn_phyloseq.rds")

sample_info <- sample_data(crohn_phyobj)
otu_mat <- otu_table(crohn_phyobj)


# ancombc
ancombc_result <- ancombc(crohn_phyobj, formula = 'Crohn', p_adj_method='BH',
                          group='Crohn', struc_zero=FALSE)

decision <- ancombc_result$res$diff_abn


ancombc_summary <- data.frame(Taxa_Name = row.names(decision),
                              Pval = ancombc_result$res$q_val$Crohn,
                              DA = ancombc_result$res$diff_abn$Crohn)


write.csv(ancombc_summary, "real_data/Crohn_Vandeputte/Crohn_ANCOMBC_result.csv",
          row.names = FALSE)



