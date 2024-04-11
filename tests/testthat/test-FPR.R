

test_that("null case", {
  null_output <- adapt(input_data=phyobj_null, 
                       cond.var="main", adj.var="confounder")
  null_results <- summary(null_output)
  FPR <- mean(null_results$pval < 0.05)
  expect_lt(FPR, 0.05)
})
