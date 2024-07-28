


#' @useDynLib ADAPT
#' @importFrom stats pchisq
#' @importFrom stats smooth.spline
#' @importFrom stats predict
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @importFrom RcppParallel RcppParallelLibs
count_ratio <- function(count_table, design_matrix, reftaxa=NULL, censor=1, test_all=FALSE){
  alltaxa <- colnames(count_table)
  if (is.null(reftaxa)) reftaxa <- alltaxa
  TBDtaxa <- NULL
  if (test_all){ # examine count ratio between other taxa and the subset
    TBDtaxa <- alltaxa
  } else{ # relative abundance within the subset of taxa
    TBDtaxa <- reftaxa
  }

  refcounts <- rowSums(count_table[, reftaxa]) # sum of reference taxa as denominator
  # remove the samples with zero counts in the reference set
  null_filter <- refcounts != 0
  refcounts <- refcounts[null_filter] # denominator
  subset_designmatrix <- design_matrix[null_filter, ]
  TBD_counts <- count_table[null_filter, TBDtaxa, drop=FALSE]
  existence <- 1*(TBD_counts > 0) # indicator matrices of zero counts
  TBD_counts[TBD_counts == 0] <- censor
  prevalences <- colMeans(existence)
  # metadata <- metadata[null_filter, ,drop=FALSE]
  # design_mat <- cbind(1, metadata[, c(covar, adjust)]) |> as.matrix()


  CR_result <- data.frame(Taxa = TBDtaxa,
                          prevalence=prevalences,
                          log10foldchange=0,
                               teststat=0,
                               pval=0)

  estimation_result <- cr_estim(count_mat=TBD_counts, refcounts=refcounts, 
                                Delta=existence, X=subset_designmatrix)
  CR_result$log10foldchange <- estimation_result[, 1] / log(10)
  CR_result$teststat <- estimation_result[, 2]
  CR_result$log10foldchange[estimation_result[, 3] == 1] <- NA # some taxa may be too rare for statistical inference
  CR_result$teststat[estimation_result[, 3] == 1] <- NA # some taxa may be too rare for statistical inference
  CR_result$pval <- 1 - pchisq(CR_result$teststat, 1)

  return(CR_result)
}
