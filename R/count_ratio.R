

#' Censored Regression for Count Ratios
#' @useDynLib PTDA
#' @param otu_table microbiome abundance matrix. All entries are integers
#' @param metadata sample metadata containing the covariate of interest and confounding factors for adjustment
#' @param covar the name of the covariate of interest
#' @param adjust the names of the variables for adjustment, NULL if no adjustments
#' @param refgenes the list of reference genes
#' @param complement If true, model the count ratios between each taxa beyond the reference set and the sum of taxa in the reference set. If false, model the relative abundance of taxa within the reference set.
#' @param boot whether to use bootstrap to estimate scaling factors for the test statistics, default TRUE
#' @param boot_replicate number of bootstrap replicates for each gene to estimate scaling factors
#' @param n_boot_gene number of genes to apply bootstrap
#' @importFrom stats pchisq
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppParallel
#' @returns The effect sizes, standard errors and p values for each taxon
#' @export
count_ratio <- function(otu_table, metadata, covar, adjust=NULL, refgenes=NULL, complement=FALSE,
                        boot=TRUE, boot_replicate=500, n_boot_gene=100){
  allgenes <- colnames(otu_table)
  if (is.null(refgenes)) refgenes <- allgenes
  TBDgenes <- NULL
  if (complement){ # examine count ratio between other taxa and the subset
    TBDgenes <- setdiff(allgenes, refgenes)
  } else{ # relative abundance within the subset of taxa
    TBDgenes <- refgenes
  }

  depths <- rowSums(otu_table[, refgenes])
  # remove the samples with zero counts in the reference set
  null_filter <- depths != 0
  depths <- depths[null_filter] # denominator
  TBD_counts <- otu_table[null_filter, TBDgenes, drop=FALSE]
  existence <- 1*(TBD_counts > 0) # indicator matrices of zero counts
  TBD_counts[!existence] <- 1
  CR_mat <- log(TBD_counts / depths)
  metadata <- metadata[null_filter, ,drop=FALSE]
  design_mat <- cbind(1, metadata[, c(covar, adjust)]) |> as.matrix()


  CR_result <- data.frame(Gene = TBDgenes,
                               effect=0,
                               teststat=0,
                               pval=0)

  estimation_result <- cr_estim(Y=CR_mat, Delta=existence, X=design_mat)
  CR_result$effect <- estimation_result[, 1]
  CR_result$teststat <- estimation_result[, 2]
  CR_result$effect[estimation_result[, 3] == 1] <- NA # some rare taxa may be difficult to estimate
  CR_result$teststat[estimation_result[, 3] == 1] <- NA # some rare taxa may be difficult to estimate

  if(boot){ # use bootstrap to estimate scaling factor
    chisq_estim <- boot_estim(Y=CR_mat, Delta=existence, X=design_mat,
                              boot_replicate=boot_replicate, n_boot_gene=n_boot_gene)
    scaling_factor <- mean(chisq_estim[, 2], na.rm=T)
    CR_result$teststat <- CR_result$teststat / scaling_factor
  }

  CR_result$pval <- 1 - pchisq(CR_result$teststat, 1)

  return(CR_result)
}
