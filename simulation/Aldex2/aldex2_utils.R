library(pROC)

# evaluate performance
evaluation <- function(taxa_truth, aldex_result, nullcase=FALSE){
  result <- NULL
  if (nullcase){
    FPR <- mean(aldex_result$wi.ep < 0.05)
    result <- list(FPR=FPR)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[rownames(aldex_result), ] == 0
    roc_obj <- roc(DA_truth_binary, aldex_result$wi.ep)
    
    aldex_DA_taxa <- rownames(aldex_result)[aldex_result$wi.eBH < 0.05]
    check_DAtaxa <- aldex_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}
