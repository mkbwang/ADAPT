

load("R/sysdata.rda")

count_mat <- simulated_null_data$otu_table
count_mat <- t(count_mat)
metadata <- simulated_null_data$metadata


depths <- rowSums(count_mat)
# remove the samples with zero counts in the reference set
null_filter <- depths != 0
depths <- depths[null_filter] # denominator
TBD_counts <- count_mat[null_filter, colnames(count_mat), drop=F]
existence <- 1*(TBD_counts > 0) # indicator matrices of zero counts
TBD_counts[!existence] <- 1
CR_mat <- log(TBD_counts / depths)
metadata <- metadata[null_filter, ,drop=FALSE]
design_mat <- cbind(1, metadata[, c("X1", "X2")]) |> as.matrix()



result <- count_ratio(otu_table=count_mat, metadata=metadata, covar="X1", adjust="X2", boot=F)

example_ID <- which.min(result$teststat)
example_ID <- 16

example <- LRT_censored_regression(Y=CR_mat[, example_ID], Delta=existence[, example_ID], X=design_mat,
                                   Firth=F, dist="lognormal")

#
# failed_genes <- result$Gene[is.na(result$effect)]


# rm(list=ls())
# testcase <- testcasesim(ngene=90, nsample=50, prev=0.4, seed=2030)
#
# ptm <- proc.time()
# real_estimates <- cr_estim(testcase$Ymat, testcase$Delta, testcase$Xmat)
# proc.time() - ptm
#
# real_estimates <- cbind(seq(1, nrow(real_estimates)),
#                         real_estimates)
# colnames(real_estimates) <- c("ID", "Estimate", "Teststat", "Fail")
#
#
# ptm <- proc.time()
# boot_chisqs <- boot_estim(testcase$Ymat, testcase$Delta, testcase$Xmat,
#                           boot_replicate = 500)
# proc.time() - ptm
#
# colnames(boot_chisqs) <- c("ID", "Chisq")

# result_mat <- outcome$inference_result
# sum(result_mat[, 3])
#
# failed_indices <- which(result_mat[, 3] == 1)
#
# prevalences <- colMeans(large_testcase$Delta[, failed_indices])
#
# Y_example <- large_testcase$Ymat[, 4]
# delta_example <- large_testcase$Delta[, 4]
# X_example <- large_testcase$Xmat
#
combined_data <- data.frame(CR_mat[, example_ID], existence[, example_ID], design_mat[, -1])
sum(combined_data[, 2])

write.table(combined_data, file=file.path("/home/wangmk/UM/Research/MDAWG/Bartlett-LRT", "example_data.txt"), col.names=F,
            row.names=F)
