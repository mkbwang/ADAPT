
source('POLDA/overdisperse_GLM.R')
source('POLDA/Reference_Selection.R')

library(ACAT)
cauchy_combination <- function(pvalvec){
  pvalvec <- na.omit(pvalvec)
  return(ACAT(pvalvec))
}


polda <- function(otu_table, metadata, covar,
                  covartype=c("categorical", "numerical"),
                  startdrop=c("median", "mean")){

  covartype <- match.arg(covartype)
  taxa_names <- row.names(otu_table) # assume taxa are rows
  struct_zero_taxa <- c()
  if (covartype == "categorical"){
    # get rid of taxa with structural zeros
    uniq_categories <- unique(metadata[, covar])
    numcategories <- length(uniq_categories)
    otu_iszero <- otu_table == 0
    otu_structzero <- matrix(FALSE, nrow=nrow(otu_iszero), ncol=numcategories)
    for (j in 1:numcategories){
      otu_iszero_subset <- otu_iszero[, metadata[, covar] == uniq_categories[j]]
      otu_structzero[, j] <- apply(otu_iszero_subset, 1, all)
    }
    is_struct_zero <- apply(otu_structzero, 1, any)
    struct_zero_taxa <- taxa_names[is_struct_zero]
    otu_table <- otu_table[!is_struct_zero, ]
    taxa_names <- taxa_names[!is_struct_zero]
  }

  pairglm_result <- pairwise_GLM(otu_table, metadata, covar) # overdispersed GLM
  reftaxa <- backward_selection(otu_table, pairglm_result, start=startdrop) # reference taxa selection
  refglm_result <-reference_GLM(otu_table, metadata, covar, reftaxa) # combine counts of reference taxa


  DiffTaxa <- refglm_result$Taxon[refglm_result$pval_adjust < 0.05]

  result <- list(Structural_zero_Taxa = struct_zero_taxa,
                 Reference_Taxa = reftaxa,
                 DA_taxa = DiffTaxa,
                 P_Value = refglm_result)

  return(result)
}

