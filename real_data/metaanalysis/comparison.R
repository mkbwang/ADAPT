
schubert_taxonomy_df <- read.csv('real_data/metaanalysis/schubert_performance.csv',
                                 row.names=1)

schubert_subset <- schubert_taxonomy_df[, c(6, 8, 9, 10)]
colnames(schubert_subset) <- c("Genus", "polda_DA_schubert", "ancombc_DA_schubert",
                               "aldex_DA_schubert")


singh_taxonomy_df <-  read.csv('real_data/metaanalysis/edd_singh_performance.csv',
                               row.names=1)
singh_subset <- singh_taxonomy_df[, c(6,8,9,10)]
colnames(singh_subset) <- c("Genus", "polda_DA_singh", "ancombc_DA_singh",
                               "aldex_DA_singh")

comparison <- schubert_subset %>% inner_join(singh_subset, by="Genus")

