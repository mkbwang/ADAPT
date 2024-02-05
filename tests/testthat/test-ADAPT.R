

test_that("null case", {
  result <- adapt(input_data=null_data, cond.var="X1", adj.var=c("X2"))
  pval_df <- result@details
  FPR <- mean(pval_df$pval < 0.05)
  expect_lt(FPR, 0.05)
})


test_that("real data", {
  real_result <-adapt(input_data=ecc16s_12, cond.var="CaseEver", adj.var=c("geo_loc_name", "host_sex"),
                 ref.cond="Control")
})

