
rm(list=ls())
# library(PTDA)
# load("R/sysdata.rda")

folder <- "experiment/sim_example_enrich"
seed <- 27
simulated_data <- readRDS(file.path(folder, sprintf("seed_%d", seed),
                                    sprintf("sim_n50_p10_d100000_f3.0_s%d.rds", seed)))



count_mat <- simulated_data$count_mat
prevalences <- rowMeans(count_mat > 0)
count_mat <- count_mat[prevalences > 0.1, ]
metadata <- simulated_data$sample_metadata
taxa_info <- simulated_data$taxa_info
true_DAtaxa <- rownames(taxa_info)[taxa_info$Logfold != 0]

result1_noboot <- ptda(otu_table = count_mat, metadata=metadata,
               covar="Group", genes_are_rows = TRUE, boot=F, prevalence_cutoff = 0.1)
pval_df_noboot <- result1_noboot$P_Value

ptm <- proc.time()
result1_boot <- ptda(otu_table = count_mat, metadata=metadata,
                     covar="Group", genes_are_rows = TRUE, boot=T, prevalence_cutoff = 0.1,
                     n_boot_gene=100)
proc.time() - ptm

pval_df_boot <- result1_boot$P_Value
count_mat <- t(count_mat)
FP_noboot <- result1_noboot$DA_Gene[which(!(result1_noboot$DA_Gene %in% true_DAtaxa))]
FP_boot <- result1_boot$DA_Gene[which(!(result1_boot$DA_Gene %in% true_DAtaxa))]

extrataxa <- setdiff(result1_boot$DA_Gene, result1_noboot$DA_Gene)
counts_FP <- count_mat[, FP_noboot]
colMeans(counts_FP > 0)
counts_extra <- count_mat[, extrataxa]
library(dplyr)
View(pval_df_noboot %>% filter(Gene %in% FP_noboot))


final_result <- count_ratio(otu_table = count_mat, metadata=metadata, covar="Group",
                            refgenes=result1_boot$Reference_Gene,
                            test_all=TRUE, boot=TRUE, boot_replicate=1000)


boot_result <- final_result$boot_estim
plot(boot_result$log_prev, boot_result$log_teststat)
abline(h=0, col="red", lwd=3, lty=2)

lm_model <- lm(log_teststat ~ log_prev, data=boot_result)
