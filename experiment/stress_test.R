
rm(list=ls())
testcase <- testcasesim(ngene=90, nsample=50, prev=0.4, seed=2030)

ptm <- proc.time()
real_estimates <- cr_estim(testcase$Ymat, testcase$Delta, testcase$Xmat)
proc.time() - ptm

real_estimates <- cbind(seq(1, nrow(real_estimates)),
                        real_estimates)
colnames(real_estimates) <- c("ID", "Estimate", "Teststat", "Fail")


ptm <- proc.time()
boot_chisqs <- boot_estim(testcase$Ymat, testcase$Delta, testcase$Xmat,
                          boot_replicate = 500)
proc.time() - ptm

colnames(boot_chisqs) <- c("ID", "Chisq")

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
# combined_data <- data.frame(Y_example, delta_example, X_example[, -1])
# sum(combined_data[, 2])
#
# write.table(combined_data, file=file.path("/home/wangmk/UM/Research/MDAWG/Bartlett-LRT", "example_data.txt"), col.names=F,
#             row.names=F)
