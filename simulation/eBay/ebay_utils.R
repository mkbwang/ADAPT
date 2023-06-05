library(pROC)


evaluation <- function(taxa_truth, ebay_result, nullcase=FALSE){
  raw_pvals <- ebay_result$final.p
  taxa_names <- names(raw_pvals)
  adjusted_pvals <- p.adjust(raw_pvals, method="BH")
  performance <- NULL
  if (nullcase){
    FPR <- mean(raw_pvals < 0.05)
    performance <- list(FPR=FPR)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[taxa_names, ] == 0
    roc_obj <- roc(DA_truth_binary, raw_pvals)
    ebay_DAtaxa <- taxa_names[adjusted_pvals < 0.05]
    check_DAtaxa <- ebay_DAtaxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    performance <- list(FDR = FDR,
                        Power=Power,
                        AUROC = roc_obj$auc)
  }
  result <- list(performance=performance, pvals=raw_pvals, adjusted_pvals=adjusted_pvals)
  return(result)
}

