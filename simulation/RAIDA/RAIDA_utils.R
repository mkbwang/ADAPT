library(pROC)

# reorganize columns so that samples from the same group are next to each other
reorganize <- function(AGP_data){
  countmat <- AGP_data$count_mat
  metadata <- AGP_data$sample_metadata
  countmat_1 <- countmat[, metadata$X == 0]
  countmat_2 <- countmat[, metadata$X == 1]
  reorganized_countmat <- cbind(countmat_1, countmat_2)
  return(data.frame(reorganized_countmat))
}

# evaluation
evaluation <- function(taxa_truth, raida_result, nullcase=FALSE){
  raida_df <- raida_result$result
  raida_pvals <- raida_df$p
  result <- NULL
  if (nullcase){
    result <- list(FPR = mean(raida_pvals < 0.05))
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[rownames(raida_df), ] == 0
    roc_obj <- roc(DA_truth_binary, raida_df$p)
    check_reftaxa <- raida_result$reference.features %in% true_DA_taxa
    all_DA_taxa <- rownames(raida_df)[raida_df$p.adj < 0.05]
    check_DAtaxa <- all_DA_taxa %in% true_DA_taxa
    reftaxa_error <- mean(check_reftaxa)
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(reftaxa_error=reftaxa_error,
                   FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  return(result)
}

