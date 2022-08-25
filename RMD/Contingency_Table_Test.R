
# load data
rm(list=ls())
setwd("/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/")
source('RMD/preprocess.R')
source('RMD/ANCOM.R')
folder = 'RMD/CAARS_data'
load(file.path(folder, 'CAARS_data_order_level.Rdata'))
processed_data <- preprocess(phylodata_order, 'SAMPLE_ID', "asthma")


filtered_count <- processed_data$feature_table
filtered_metadata <- processed_data$meta_data
struc_zero <- processed_data$structure_zeros


included_genus <- names(which(rowSums(struc_zero) == 0))
taxa_pairs <- combn(included_genus, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("Taxa1", "Taxa2")
taxa_pairs$pval <- NA

asthma_counts <- filtered_count[, filtered_metadata$asthma == 1]
nonasthma_counts <- filtered_count[, filtered_metadata$asthma == 0]


for (j in 1:nrow(taxa_pairs)){

  selected_taxon1 <- taxa_pairs[j, 1]
  selected_taxon2 <- taxa_pairs[j, 2]

  asthma_tx1 <- asthma_counts[selected_taxon1, ] %>% as.matrix()
  asthma_tx2 <- asthma_counts[selected_taxon2, ] %>% as.matrix()
  nonasthma_tx1 <- nonasthma_counts[selected_taxon1, ] %>% as.matrix()
  nonasthma_tx2 <- nonasthma_counts[selected_taxon2, ] %>% as.matrix()


  asthma_dbpos <- sum((asthma_tx1 > 0) & (asthma_tx2 > 0), na.rm=TRUE)
  asthma_dbneg <- sum((asthma_tx1 == 0) & (asthma_tx2 == 0), na.rm=TRUE)
  asthma_t1post2neg <- sum((asthma_tx1 > 0) & (asthma_tx2 == 0), na.rm=TRUE)
  asthma_t1negt2pos <- sum((asthma_tx1 == 0) & (asthma_tx2 > 0), na.rm=TRUE)

  control_dbpos <- sum((nonasthma_tx1 > 0) & (nonasthma_tx2 > 0), na.rm=TRUE)
  control_dbneg <- sum((nonasthma_tx1 == 0) & (nonasthma_tx2 == 0), na.rm=TRUE)
  control_t1post2neg <- sum((nonasthma_tx1 > 0) & (nonasthma_tx2 == 0), na.rm=TRUE)
  control_t1negt2pos <- sum((nonasthma_tx1 == 0) & (nonasthma_tx2 > 0), na.rm=TRUE)

  contingency_table <- matrix(c(asthma_dbpos, asthma_dbneg, asthma_t1post2neg, asthma_t1negt2pos,
                                control_dbpos, control_dbneg, control_t1post2neg, control_t1negt2pos),
                              nrow=4)

  testresult <- fisher.test(contingency_table)
  taxa_pairs[j, 3] <- testresult$p.value

}

taxa_pairs$adjustedp <- p.adjust(taxa_pairs$pval, method='BH')
write.csv(taxa_pairs, 'RMD/CAARS_Model_Summary/fisher_test.csv', row.names = FALSE)

