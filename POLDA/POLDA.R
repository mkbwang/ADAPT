
source('/home/wangmk/UM/Research/MDAWG/POLDA/POLDA/overdisperse_GLM.R')
source('/home/wangmk/UM/Research/MDAWG/POLDA/POLDA/Reference_Selection.R')


polda <- function(otu_table, metadata, covar,
                  covartype=c("categorical", "numerical")){

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

  cat("Estimate relative abundance ratio for all taxa...\n")
  alltaxa_result <- relabd_GLM(otu_table, metadata, covar)
  
  med_relabd_odds <- median(alltaxa_result$effect)
  cat("Backward reference taxa selection...\n")
  # backward selection for the "depleted" taxa group
  Taxa_D <- alltaxa_result$Taxa[alltaxa_result$effect < med_relabd_odds]
  while(1) {
    otu_table_D <- otu_table[Taxa_D, ]
    deplete_taxa_result <- relabd_GLM(otu_table_D, metadata, covar)
    if (all(deplete_taxa_result$teststat > qnorm(0.05)) & 
        all(deplete_taxa_result$teststat < qnorm(0.95))) break
    Taxa_D <- deplete_taxa_result$Taxa[deplete_taxa_result$effect > 0]
  }
  
  # backward selection for the "enriched" taxa group
  Taxa_E <- alltaxa_result$Taxa[alltaxa_result$effect > med_relabd_odds]
  while(1){
    otu_table_E <- otu_table[Taxa_E, ]
    enriched_taxa_result <- relabd_GLM(otu_table_E, metadata, covar)
    if (all(enriched_taxa_result$teststat > qnorm(0.05)) & 
        all(enriched_taxa_result$teststat < qnorm(0.95))) break
    Taxa_E <- enriched_taxa_result$Taxa[enriched_taxa_result$effect < 0]
  }
  
  reftaxa <- c(Taxa_D, Taxa_E)
  
  # pairglm_result <- pairwise_GLM(otu_table, metadata, covar) # overdispersed GLM
  # selection_result <- backward_selection(otu_table, pairglm_result, start=startdrop) # reference taxa selection
  # reftaxa <- selection_result$reftaxa
  cat("Find differentially abundant taxa based on reference taxa...\n")
  refglm_result <-reference_GLM(otu_table, metadata, covar, reftaxa) # combine counts of reference taxa


  DiffTaxa <- refglm_result$Taxon[refglm_result$pval_adjust < 0.05]

  result <- list(Structural_zero_Taxa = struct_zero_taxa,
                 Reference_Taxa = reftaxa,
                 DA_taxa = DiffTaxa,
                 P_Value = refglm_result)

  return(result)
}

