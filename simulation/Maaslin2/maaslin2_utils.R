library(pROC)

evaluation <- function(taxa_truth, maaslin2_result, nullcase=FALSE){
  maaslin2_core <- maaslin2_result$results[, c("feature", "pval", "qval")]
  raw_pvals <- maaslin2_core$pval
  adjusted_pvals <- maaslin2_core$qval
  taxa_names <- maaslin2_core$feature
  performance <- NULL
  if (nullcase){
    FPR <- mean(raw_pvals < 0.05)
    performance <- list(FPR=FPR)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[taxa_names, ] == 0
    roc_obj <- roc(DA_truth_binary, raw_pvals)
    maaslin2_DAtaxa <- taxa_names[adjusted_pvals < 0.05]
    check_DAtaxa <- maaslin2_DAtaxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    performance <- list(FDR = FDR,
                        Power=Power,
                        AUROC = roc_obj$auc)
  }
  result <- list(performance=performance, pvals=maaslin2_core)
  return(result)
}
