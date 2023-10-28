

test_that("null case", {
  count_mat <- simulated_null_data$otu_table
  metadata <- simulated_null_data$metadata
  result2_boot <- adapt(otu_table = count_mat, metadata=metadata,
                       covar="X1", adjust="X2", genes_are_rows = TRUE, boot=T)
  pval_df <- result2_boot$P_Value
  FPR <- mean(pval_df$pval < 0.05)
  expect_lt(FPR, 0.05)
})


