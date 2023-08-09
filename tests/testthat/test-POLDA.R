test_that("null case", {
  polda_null <- polda(otu_table=example_null$count_mat,
                      metadata=example_null$sample_metadata,
                      covar="X", ratio_model="loglogistic")
  polda_pvals <- polda_null$P_Value
  expect_lt(mean(polda_pvals$pval < 0.05), 0.05)
})


test_that("DA case", {
  polda_DA <- polda(otu_table=example_DA$count_mat,
                    metadata=example_DA$sample_metadata,
                    covar="X", ratio_model="loglogistic")
  taxa_truth <- example_DA$taxa_info
  true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
  check_DAtaxa <- polda_DA$DA_taxa %in% true_DA_taxa
  check_reftaxa <- polda_DA$Reference_Taxa %in% true_DA_taxa

  reftaxa_error <- mean(check_reftaxa)
  expect_equal(reftaxa_error, 0)
  FDR <- 1 - mean(check_DAtaxa)
  expect_lt(FDR, 0.05)
  Power <- sum(check_DAtaxa) / length(true_DA_taxa)
  expect_gt(Power, 0.5)
})
