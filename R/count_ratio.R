

#' Censored Regression for Count Ratios
#' @useDynLib ADAPT
#' @param otu_table microbiome abundance matrix. All entries are integers
#' @param metadata sample metadata containing the covariate of interest and confounding factors for adjustment
#' @param covar the name of the covariate of interest
#' @param adjust the names of the variables for adjustment, NULL if no adjustments
#' @param reftaxa the list of reference taxa
#' @param test_all If true, model the count ratios between all taxa against the reference set; otherwise, only model the relative abundance within the reference set
#' @param boot whether to use bootstrap to estimate scaling factors for the test statistics, default TRUE
#' @param boot_replicate number of bootstrap replicates for each gene to estimate scaling factors
#' @param n_boot_taxa number of taxa to apply bootstrap
#' @importFrom stats pchisq
#' @importFrom stats smooth.spline
#' @importFrom stats predict
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppParallel
#' @returns The effect sizes, standard errors and p values for each taxon
#' @export
count_ratio <- function(otu_table, metadata, covar, adjust=NULL, reftaxa=NULL, test_all=FALSE,
                        boot=TRUE, boot_replicate=2000, n_boot_taxa=500){
  alltaxa <- colnames(otu_table)
  if (is.null(reftaxa)) reftaxa <- alltaxa
  TBDtaxa <- NULL
  if (test_all){ # examine count ratio between other taxa and the subset
    TBDtaxa <- alltaxa
  } else{ # relative abundance within the subset of taxa
    TBDtaxa <- reftaxa
  }

  refcounts <- rowSums(otu_table[, reftaxa]) # sum of reference taxa as denominator
  # remove the samples with zero counts in the reference set
  null_filter <- refcounts != 0
  refcounts <- refcounts[null_filter] # denominator
  TBD_counts <- otu_table[null_filter, TBDtaxa, drop=FALSE]
  existence <- 1*(TBD_counts > 0) # indicator matrices of zero counts
  prevalences <- colMeans(existence)
  metadata <- metadata[null_filter, ,drop=FALSE]
  design_mat <- cbind(1, metadata[, c(covar, adjust)]) |> as.matrix()


  CR_result <- data.frame(Taxa = TBDtaxa,
                          prevalence=prevalences,
                               effect=0,
                               teststat=0,
                               pval=0)

  estimation_result <- cr_estim(count_mat=TBD_counts, refcounts=refcounts, 
                                Delta=existence, X=design_mat)
  CR_result$effect <- estimation_result[, 1]
  CR_result$teststat <- estimation_result[, 2]
  CR_result$effect[estimation_result[, 3] == 1] <- NA # some taxa may be too rare for statistical inference
  CR_result$teststat[estimation_result[, 3] == 1] <- NA # some taxa may be too rare for statistical inference

  chisq_estim_df <- NULL
  if(boot){ # use bootstrap to estimate scaling factor
    chisq_estim <- boot_estim(count_mat=TBD_counts, refcounts=refcounts, Delta=existence, X=design_mat,
                              boot_replicate=boot_replicate, n_boot_taxa=n_boot_taxa)
    
    scaling_factor <- chisq_estim[, 2]
    if (n_boot_taxa < length(alltaxa)){ # chosen taxa to bootstrap is a subset of all taxa

      selected_taxa <- chisq_estim[, 1] # column numbers of selected taxa for scale estimation
      selected_prevalences <- colMeans(existence[, selected_taxa])
      
      # fit a cubic spline to impute the scaling factors for the unselected taxa based on prevalence
      ssfit <- smooth.spline(selected_prevalences, log(scaling_factor))
      predict_ssfit <- predict(ssfit, CR_result$prevalence)
      scaling_factor <- exp(predict_ssfit$y)
      scaling_factor[selected_taxa] <- chisq_estim[, 2] # don't impute those that have been estimated
      
    }
    CR_result$teststat <- CR_result$teststat / scaling_factor

  }

  CR_result$pval <- 1 - pchisq(CR_result$teststat, 1)

  return(list(CR_result=CR_result,
              boot_estim=chisq_estim_df))
}
