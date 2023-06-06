library(pROC)
library(GUniFrac)
library(DESeq2)

# fit model
gmpr_deseq <- function(countmat, metadata){
  GMPR_factors <- GMPR(OTUmatrix=countmat)
  dds<- DESeq2::DESeqDataSetFromMatrix(countData=countmat, 
                                            colData=metadata, 
                                            design=stats::as.formula("~X"))
  
  sizeFactors(dds) <- GMPR_factors
  deseq_out <- DESeq(object=dds, test = "LRT", reduced = ~1)
  res_df <- results(deseq_out) |> as.data.frame()
  return(res_df)
}


# evaluate model
evaluation <- function(taxa_truth, deseq_result, nullcase=F){
  
  raw_pvals <- deseq_result$pvalue
  adjusted_pvals <- deseq_result$padj
  taxa_names <- rownames(deseq_result)
  performance <- NULL
  
  if (nullcase){
    FPR <- mean(raw_pvals < 0.05)
    performance <- list(FPR=FPR)
  } else{
    true_DA_taxa <- rownames(taxa_truth)[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[taxa_names, ] == 0
    roc_obj <- roc(DA_truth_binary, raw_pvals)
    deseq_DAtaxa <- taxa_names[adjusted_pvals < 0.05 & !is.na(adjusted_pvals)]
    check_DAtaxa <- deseq_DAtaxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    performance <- list(FDR = FDR,
                        Power=Power,
                        AUROC = roc_obj$auc)
  }
  output <- list(performance=performance, pval=deseq_result)
  return(output)
}
