
rm(list=ls())
load("R/sysdata.rda")

count_mat <- simulated_null_data$otu_table
metadata <- simulated_null_data$metadata

result1_noboot <- ptda(otu_table = count_mat, metadata=metadata,
               covar="X1", genes_are_rows = TRUE, boot=F)

result1_boot <- ptda(otu_table = count_mat, metadata=metadata,
                     covar="X1", genes_are_rows = TRUE, boot=T)

result2_noboot <- ptda(otu_table = count_mat, metadata=metadata,
                       covar="X1", adjust="X2", genes_are_rows = TRUE, boot=F)

ptm <- proc.time()
result2_boot <- ptda(otu_table = count_mat, metadata=metadata,
                     covar="X1", adjust="X2", genes_are_rows = TRUE, boot=T)
proc.time() - ptm
