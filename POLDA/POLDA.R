
source('/home/wangmk/MDAWG/POLDA/POLDA/overdisperse_GLM.R')

library(ClassComparison) # for fitting BUM to p values

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

  
  reftaxa <- taxa_names # initially all the taxa are reference taxa(relative abundance)
  while(1){
   
    relabd_result <- reference_GLM(count_data = otu_table, metadata = metadata, 
                                    covar, reftaxa = reftaxa, complement=FALSE)
    
    estimated_effect <- relabd_result$effect
    pvals <- relabd_result$pval
    names(pvals) <- relabd_result$Taxon
    names(estimated_effect) <- relabd_result$Taxon
    # check distribution of p values
    bumfit <- Bum(pvals)
    lambda_hat <- bumfit@lhat
    a_hat <- bumfit@ahat
    # check if the sum of log likelihood is larger than 4 (AIC checking)
    loglik <- sum(log(likelihoodBum(bumfit)))
    if (loglik > 2){ # need to continue shrinking reference taxa set
      distance2med <- abs(estimated_effect - median(estimated_effect))
      sorted_distance <- sort(distance2med)
      ordered_taxanames <- names(sorted_distance)
      reftaxa <- ordered_taxanames[1:(length(ordered_taxanames)/2)]
    } else{
      break
    }
  }
  
  if (length(reftaxa) < length(taxa_names)){
    # Fit overdispersed GLM for all the other taxa outside reference set
    complement_result <- reference_GLM(count_data = otu_table, metadata = metadata, 
                                       covar, reftaxa = reftaxa, complement=TRUE)
    # combine p values for all the taxa
    all_GLM_results <- rbind(relabd_result, complement_result)
  } else{ # relative abundance is good enough for DAA
    all_GLM_results <- relabd_result
  }
  
  all_pvals <- all_GLM_results$pval
  names(all_pvals) <- all_GLM_results$Taxon
  # find p value cutoff for 0.05 FDR
  bumfit <- Bum(all_pvals)
  p_cutoff <- cutoffSignificant(bumfit, alpha=0.05, by="FDR")
  significant_pvals <- all_pvals[all_pvals < p_cutoff] 
  if (length(significant_pvals) > 0){
    DiffTaxa <- names(significant_pvals)
  } else{
    DiffTaxa <- c()
  }
  
  result <- list(Structural_zero_Taxa = struct_zero_taxa,
                 Reference_Taxa = reftaxa,
                 DA_taxa = DiffTaxa,
                 P_Value = all_GLM_results,
                 Pval_Cutoff = p_cutoff)

  return(result)
}

