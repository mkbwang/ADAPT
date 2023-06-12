library(pROC)
library(metagenomeSeq)

# set up MRexperiment object
setExperiment <- function(metadata, countmat){
  phenotype <- AnnotatedDataFrame(metadata)
  MR_object <- newMRexperiment(counts=countmat, phenoData=phenotype, featureData=NULL)
  MR_object <- wrenchNorm(MR_object, condition=MR_object$X)
  return(MR_object)
}

# evaluate performance
evaluation <- function(taxa_truth, MR_result, nullcase=FALSE){
  performance <- NULL
  MR_coefficients <- MRfulltable(MR_result, number=nrow(taxa_truth), adjustMethod="BH")
  MR_coefficients <- MR_coefficients[, c("logFC", "se", "pvalues", "adjPvalues")]
  if (nullcase){
    FPR <- mean(MR_coefficients$pvalues < 0.05,  na.rm=TRUE)
    performance <- list(FPR=FPR)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[rownames(MR_coefficients), ] == 0
    roc_obj <- roc(DA_truth_binary, MR_coefficients$pvalues, na.rm=TRUE)
    
    MR_DA_taxa <- rownames(MR_coefficients)[MR_coefficients$adjPvalues < 0.05]
    check_DAtaxa <- MR_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    performance <- list(FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  model_summary <- list(performance=performance, coefficients=MR_coefficients)
  return(model_summary)
}
