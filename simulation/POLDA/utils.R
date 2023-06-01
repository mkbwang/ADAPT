library(pROC)

# evaluate performance
evaluation <- function(taxa_truth, polda_result, nullcase=FALSE){
  polda_pvals <- polda_result$P_Value
  result <- NULL
  if (nullcase){
    result <- list(FPR = mean(polda_pvals$pval < 0.05))
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[polda_pvals$Taxon, ] == 0
    roc_obj <- roc(DA_truth_binary, polda_pvals$pval)
    check_reftaxa <- polda_result$Reference_Taxa %in% true_DA_taxa
    all_DA_taxa <- c(polda_result$Structural_zero_Taxa, polda_result$DA_taxa)
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


# likelihood of BUM
# llk_func <- function(lambda_hat, a_hat, p_val){
#   result <- lambda_hat + (1-lambda_hat) * a_hat * (p_val)^(a_hat - 1)
#   return(result)
# }
# Cutoff calculation for FDR
# cutoff_func <- function(lambda_hat, a_hat, alpha){
#   pi_hat <- lambda_hat + (1 - lambda_hat) * a_hat
#   result <- ((pi_hat - alpha*lambda_hat)/(alpha*(1-lambda_hat)))^(1/(a_hat - 1))
#   return(result)
# }

