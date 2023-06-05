library(pROC)



# evaluate performance
evaluation <- function(taxa_truth, ancombc_result, nullcase=FALSE){
  result <- NULL
  if (nullcase){
    pval_df <- ancombc_result$res$p_val
    FPR <- mean(pval_df[, 3] < 0.05)
    result <- list(FPR=FPR)
  } else{
    pval_df <- ancombc_result$res$p_val
    qval_df <- ancombc_result$res$q_val
    decision <- ancombc_result$res$diff_abn
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[pval_df$taxon, ] == 0
    roc_obj <- roc(DA_truth_binary, as.vector(pval_df[, 3]))
    
    ancombc_DA_taxa <- decision$taxon[decision[, 3]]
    check_DAtaxa <- ancombc_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


