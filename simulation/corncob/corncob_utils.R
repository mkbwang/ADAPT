library(pROC)


evaluation <- function(taxa_truth, corncob_result, nullcase=FALSE){
  raw_pvals <- corncob_result$p
  adjusted_pvals <- corncob_result$p_fdr
  pval_df <- cbind(raw_pvals, adjusted_pvals) |> as.data.frame()
  performance <- NULL
  if (nullcase){
    FPR <- mean(raw_pvals < 0.05)
    performance <- list(FPR=FPR)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[rownames(pval_df), ] == 0
    roc_obj <- roc(DA_truth_binary, raw_pvals)
    corncob_DAtaxa <- corncob_result$significant_taxa
    check_DAtaxa <- corncob_DAtaxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    performance <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  result <- list(pvals=pval_df, performance=performance)
  return(result)
}
