library(GUniFrac)
library(edgeR)
library(pROC)


# fit model
gmpr_edger <- function(countmat, metadata){
  GMPR_factors<- GMPR(OTUmatrix = countmat)
  dge <- edgeR::DGEList(counts = countmat)
  dge$samples$norm.factors <- GMPR_factors
  design <- stats::model.matrix(stats::as.formula("~ X"), metadata)
  dge <- edgeR::estimateDisp(dge, design)
  glmFit <- edgeR::glmFit(dge, design)
  glmRes <- edgeR::glmLRT(glmFit, coef = 2)
  return(glmRes$table)
}


# model evaluation
evaluation <- function(taxa_truth, edger_result, nullcase=F){
  raw_pvals <- edger_result$PValue
  adjusted_pvals <- p.adjust(raw_pvals, "BH")
  taxa_names <- rownames(edger_result)
  performance <- NULL
  
  if (nullcase){
    FPR <- mean(raw_pvals < 0.05)
    performance <- list(FPR=FPR)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[taxa_names, ] == 0
    roc_obj <- roc(DA_truth_binary, raw_pvals)
    edger_DAtaxa <- taxa_names[adjusted_pvals < 0.05]
    check_DAtaxa <- edger_DAtaxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    performance <- list(FDR = FDR,
                        Power=Power,
                        AUROC = roc_obj$auc)
  }
  output <- list(performance=performance, pval=edger_result)
  return(output)
}

