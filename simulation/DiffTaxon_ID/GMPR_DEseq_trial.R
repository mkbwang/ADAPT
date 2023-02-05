rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
# folder = '.'
library(dplyr)
library(GUniFrac)
library(DESeq2)


ID <- 1
input_fname <- sprintf('simulated_data_%d.rds', ID)
simulated_data <- readRDS(file.path(folder, 'data', 'semiparametric', input_fname))


count_mat <- simulated_data$otu.tab.sim

metadata <- as.data.frame(simulated_data$covariate)
colnames(metadata) <- 'X'
metadata$X <- as.factor(metadata$X)

depths <- colSums(count_mat)
normalization_factor <- GMPR(count_mat)
normalized_count <- t(t(count_mat) / normalization_factor) %>% round()


dds <- DESeqDataSetFromMatrix(countData = normalized_count,
                              colData = metadata,
                              design = ~X)

dds <- DESeq(dds)
dds_result <- results(dds)

deseq_result <- dds_result$padj < 0.05
result <- data.frame(DA = deseq_result,
                     truth = simulated_data$diff.otu.ind)
row.names(result) <- row.names(count_mat)

output_fname <- sprintf('GMPR_DESEQ2_%d.csv', ID)
write.csv(result, file.path(folder, 'GMPR_DEseq_result', output_fname))

