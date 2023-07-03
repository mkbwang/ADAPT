

#' Pooling overdispersed logistic regressions for microbiome differential abundance analysis
#'
#' @param otu_table microbiome abundance table generated from 16S rRNA sequencing or shotgun metagenomic sequencing. All entries are integers
#' @param metadata sample metadata dataframe
#' @param covar the name of the covariate of interest
#' @param adjust the names of confounders to adjust for, NULL if no adjustment
#' @param prevalence_cutoff features whose prevalence are smaller than the cutoff will be excluded from analysis
#' @param depth_cutoff a sample would be discarded if its sequencing depth is smaller than the threshold
#' @param features_are_rows whether the microbiome features are on the column or the row of the abundance table
#' @param alpha the cutoff of the adjusted p values
#' @importFrom ClassComparison Bum
#' @importFrom ClassComparison likelihoodBum
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @returns the reference feature set, the identified DA features and the p values for all the features
#' @export
polda <- function(otu_table, metadata, covar, adjust=NULL,
                  prevalence_cutoff=0.1, depth_cutoff=1000, features_are_rows=TRUE,
                  alpha=0.05){

  if(!features_are_rows) otu_table <- t(otu_table)

  prevalences <- rowMeans(otu_table != 0)
  depths <- colSums(otu_table)
  otu_table_filtered <- otu_table[prevalences > prevalence_cutoff, depths > depth_cutoff]
  taxa_names <- row.names(otu_table_filtered)
  reftaxa <- taxa_names # initially all the taxa are reference taxa(relative abundance)
  while(1){
    relabd_result <- GLM_count_ratio(count_data = otu_table_filtered, metadata = metadata,
                                    covar=covar, adjust=adjust, reftaxa = reftaxa, complement=FALSE)

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
    complement_result <- GLM_count_ratio(count_data = otu_table_filtered, metadata = metadata,
                                       covar=covar, adjust=adjust, reftaxa = reftaxa, complement=TRUE)
    # combine p values for all the taxa
    all_GLM_results <- rbind(relabd_result, complement_result)
  } else{ # relative abundance is good enough for DAA
    all_GLM_results <- relabd_result
  }

  all_pvals <- all_GLM_results$pval
  names(all_pvals) <- all_GLM_results$Taxon
  # find p value cutoff for 0.05 FDR
  all_adjusted_pvals <- p.adjust(all_pvals, method="BH")
  all_GLM_results$adjusted_pval <- all_adjusted_pvals
  # bumfit <- Bum(all_pvals)
  # p_cutoff <- cutoffSignificant(bumfit, alpha=0.05, by="FDR")
  significant_pvals <- all_adjusted_pvals[all_adjusted_pvals < alpha]
  if (length(significant_pvals) > 0){
    DiffTaxa <- names(significant_pvals)
  } else{
    DiffTaxa <- c()
  }

  result <- list(Reference_Taxa = reftaxa,
                 DA_taxa = DiffTaxa,
                 P_Value = all_GLM_results)

  return(result)
}

