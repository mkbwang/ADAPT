

#' Censored Regression for Count Ratios
#' @useDynLib PTDA
#' @param otu_table microbiome abundance matrix. All entries are integers
#' @param metadata sample metadata containing the covariate of interest and confounding factors for adjustment
#' @param covar the name of the covariate of interest
#' @param adjust the names of the variables for adjustment, NULL if no adjustments
#' @param refgenes the list of reference genes
#' @param test_all If true, model the count ratios between all taxa against the reference set; otherwise, only model the relative abundance within the reference set
#' @param boot whether to use bootstrap to estimate scaling factors for the test statistics, default TRUE
#' @param boot_replicate number of bootstrap replicates for each gene to estimate scaling factors
#' @param n_boot_gene number of genes to apply bootstrap
#' @importFrom stats pchisq
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppParallel
#' @returns The effect sizes, standard errors and p values for each taxon
#' @export
count_ratio <- function(otu_table, metadata, covar, adjust=NULL, refgenes=NULL, test_all=FALSE,
                        boot=TRUE, boot_replicate=500, n_boot_gene=100){
  allgenes <- colnames(otu_table)
  if (is.null(refgenes)) refgenes <- allgenes
  TBDgenes <- NULL
  if (test_all){ # examine count ratio between other taxa and the subset
    TBDgenes <- allgenes
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
  prevalences <- colMeans(existence)
  CR_mat <- log(TBD_counts / depths)
  metadata <- metadata[null_filter, ,drop=FALSE]
  design_mat <- cbind(1, metadata[, c(covar, adjust)]) |> as.matrix()


  CR_result <- data.frame(Gene = TBDgenes,
                          prevalence=prevalences,
                               effect=0,
                               teststat=0,
                               pval=0)

  estimation_result <- cr_estim(Y=CR_mat, Delta=existence, X=design_mat)
  CR_result$effect <- estimation_result[, 1]
  CR_result$teststat <- estimation_result[, 2]
  CR_result$effect[estimation_result[, 3] == 1] <- NA # some rare taxa may be difficult to estimate
  CR_result$teststat[estimation_result[, 3] == 1] <- NA # some rare taxa may be difficult to estimate

  chisq_estim_df <- NULL
  if(boot){ # use bootstrap to estimate scaling factor
    chisq_estim <- boot_estim(Y=CR_mat, Delta=existence, X=design_mat,
                              boot_replicate=boot_replicate, n_boot_gene=n_boot_gene)
    chisq_estim <- na.omit(chisq_estim)
    scaling_factor <- mean(chisq_estim[, 2])
    # fit linear model between log teststatistic and log prevalence to decide scaling factors
    selected_genes <- chisq_estim[, 1]
    selected_prevalences <- colMeans(existence[, selected_genes])
    chisq_estim_df <- data.frame(Gene = TBDgenes[selected_genes],
                                 teststat = chisq_estim[, 2],
                                 prev = selected_prevalences)
    chisq_estim_df$log_teststat <- log(chisq_estim_df$teststat)
    chisq_estim_df$log_prev <- log(chisq_estim_df$prev)
    regression <- lm(log_teststat ~ log_prev, data=chisq_estim_df)
    lm_coefs <- summary(regression)$coefficients
    scaling_factor <- mean(chisq_estim[, 2], na.rm=T) # default average

    if (lm_coefs[2,4] < 0.05){# scaling factors need to be different based on prevalences
        scaling_factor <- exp(lm_coefs[1,1] + lm_coefs[2,1] * log(CR_result$prevalence))
    }

    CR_result$teststat <- CR_result$teststat / scaling_factor

  }

  CR_result$pval <- 1 - pchisq(CR_result$teststat, 1)

  return(list(CR_result=CR_result,
              boot_estim=chisq_estim_df))
}
