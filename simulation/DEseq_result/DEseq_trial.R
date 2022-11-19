rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
# folder = '.'
library(dplyr)

library(foreach)
library(doParallel)
cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

library(DESeq2)

deseq_summary <- foreach(fnum=1:100, .combine=rbind,
        .packages=c('dplyr', 'DESeq2'), .inorder=FALSE) %dopar% {

  input_fname <- sprintf('simulated_data_%d.rds', fnum)
  data <- readRDS(file.path(folder, 'data', input_fname))

  abn_info <- data$mean.eco.abn
  counts <- data$obs.abn

  metadata <- data.frame(ID = colnames(counts),
                         group = rep(c(0, 1), each=50))
  metadata$group <- as.factor(metadata$group)



  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~group)

  dds <- DESeq(dds)

  dds_result <- results(dds)

  truth <- abn_info$effect.size != 1
  deseq_result <- dds_result$padj < 0.05
  performance <- table(truth, deseq_result) %>% as.vector()
  return(performance)
}


deseq_perform_df <- as.data.frame(deseq_summary)
colnames(deseq_perform_df) <- c("TN", "FN", "FP", "TP")
deseq_perform_df$FDR <- deseq_perform_df$FP /
  (deseq_perform_df$FP + deseq_perform_df$TP)
deseq_perform_df$Power <- deseq_perform_df$TP/
  (deseq_perform_df$TP + deseq_perform_df$FN)

write.csv(deseq_perform_df, file.path(folder, 'DESEQ2_performance.csv'),
          row.names=FALSE)

