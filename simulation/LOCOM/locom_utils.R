library(pROC)



# evaluate performance
evaluation <- function(taxa_truth, locom_result, nullcase=FALSE){
  result <- NULL
  if (nullcase){
    raw_pval <- as.vector(locom_result$p.otu)
    FPR <- mean(raw_pval < 0.05)
    result <- list(FPR=FPR)
  } else{
    raw_pval <- locom_result$p.otu
    taxa_names <- colnames(raw_pval)
    adjusted_pval <- as.vector(locom_result$q.otu)
    all_taxa <- rownames(taxa_truth)
    locom_DA_taxa <- locom_result$detected.otu
    true_DA_taxa <- all_taxa[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[taxa_names, ] == 0
    roc_obj <- roc(DA_truth_binary, as.vector(raw_pval))
    check_DAtaxa <- locom_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


