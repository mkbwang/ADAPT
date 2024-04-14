

test_that("balanced change", {
  phyobj_balanced <- balanced_case$phyobj
  balanced_truth <- balanced_case$truth
  balanced_output <- adapt(input_data=phyobj_balanced, 
                           cond.var="main", adj.var=c("confounder"))
  balanced_results <- summary(balanced_output, 'da')
  signals <- balanced_results$Taxa
  true_signals <- balanced_truth$taxaname[balanced_truth$isDA]
  FDR_balanced <- 1-mean(signals %in% true_signals)
  expect_lt(FDR_balanced, 0.05)
})


test_that("unbalanced change", {
  phyobj_unbalanced <- unbalanced_case$phyobj
  unbalanced_truth <- unbalanced_case$truth
  unbalanced_output <- adapt(input_data=phyobj_unbalanced, 
                             cond.var="main", adj.var="confounder")
  unbalanced_result <- summary(unbalanced_output, 'da')
  signals <- unbalanced_result$Taxa
  true_signals <- unbalanced_truth$taxaname[unbalanced_truth$isDA]
  FDR_unbalanced <- 1-mean(signals %in% true_signals)
  expect_lt(FDR_unbalanced, 0.05)
})

