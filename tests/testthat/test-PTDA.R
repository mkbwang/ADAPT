

test_that("null case", {

  # metadata <- null_example$metadata
  # count_mat <- null_example$count_mat
  # ptda_null <- ptda(otu_table=count_mat,
  #                       metadata=metadata,
  #                       covar="X", ratio_model="lognormal", zero_censor=1)
  # ptda_pvals <- ptda_null$P_Value
  # FPR <- mean(ptda_pvals$pval < 0.05)

  expect_lt(0.05, 0.06)
})


# test_that("DA case 1", {
#   polda_DA <- polda(otu_table=example_DA_1$count_mat,
#                     metadata=example_DA_1$sample_metadata,
#                     covar="X", ratio_model="loglogistic")
#   taxa_truth <- example_DA_1$taxa_info
#   true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
#   check_DAtaxa <- polda_DA$DA_taxa %in% true_DA_taxa
#   check_reftaxa <- polda_DA$Reference_Taxa %in% true_DA_taxa
#
#   reftaxa_error <- mean(check_reftaxa)
#   expect_equal(reftaxa_error, 0)
#   FDR <- 1 - mean(check_DAtaxa)
#   expect_lt(FDR, 0.05)
#   Power <- sum(check_DAtaxa) / length(true_DA_taxa)
#   expect_gt(Power, 0.5)
# })
#
#
# test_that("DA case 2", {
#   polda_DA <- polda(otu_table=example_DA_2$count_mat,
#                     metadata=example_DA_2$sample_metadata,
#                     covar="X", ratio_model="loglogistic")
#   taxa_truth <- example_DA_2$taxa_info
#   true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
#   check_DAtaxa <- polda_DA$DA_taxa %in% true_DA_taxa
#   check_reftaxa <- polda_DA$Reference_Taxa %in% true_DA_taxa
#
#   reftaxa_error <- mean(check_reftaxa)
#   expect_equal(reftaxa_error, 0)
#   FDR <- 1 - mean(check_DAtaxa)
#   expect_lt(FDR, 0.05)
#   Power <- sum(check_DAtaxa) / length(true_DA_taxa)
#   expect_gt(Power, 0.5)
# })
