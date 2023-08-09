

#' Censored Regression for Count Ratios
#'
#' @param count_data microbiome abundance matrix. All entries are integers
#' @param metadata sample metadata containing the covariate of interest and confounding factors for adjustment
#' @param covar the name of the covariate of interest
#' @param adjust the names of the variables for adjustment, NULL if no adjustments
#' @param reftaxa the list of reference taxa
#' @param complement If true, model the count ratios between each taxa beyond the reference set and the sum of taxa in the reference set. If false, model the relative abundance of taxa within the reference set.
#' @param ratio_model Use loglogistic model or the Weibull model to model count ratios with excessive zeros
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom stats pnorm
#' @importFrom stats optim
#' @importFrom stats pchisq
#' @importFrom stats var
#' @returns The effect sizes, standard errors and p values for each taxon
#' @export
CR_count_ratio <- function(count_data, metadata, covar, adjust=NULL, reftaxa=NULL, complement=FALSE,
                        ratio_model=c("loglogistic", "weibull")){
  alltaxa <- rownames(count_data)
  if (is.null(reftaxa)) reftaxa <- alltaxa
  TBDtaxa <- NULL
  if (complement){ # examine count ratio between other taxa and the subset
    TBDtaxa <- setdiff(alltaxa, reftaxa)
  } else{ # relative abundance within the subset of taxa
    TBDtaxa <- reftaxa
  }

  refcounts <- colSums(count_data[reftaxa, ])
  # remove the samples with zero counts in the reference set
  null_filter <- refcounts != 0
  refcounts <- refcounts[null_filter]
  TBD_counts <- count_data[TBDtaxa, null_filter, drop=FALSE]
  metadata <- metadata[null_filter, ,drop=FALSE]
  num_taxa <- length(TBDtaxa)
  CR_result <- data.frame(Taxon = TBDtaxa,
                               effect=0,
                               teststat=0,
                               pval=0)

  cores <- availableCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)

  outcome <- foreach(j=1:nrow(CR_result), .combine=rbind,
                     .inorder=FALSE) %dopar% {

                       test_taxon <- as.character(CR_result$Taxon[j])

                       numerator <- TBD_counts[test_taxon, ] # the counts of taxa of interest
                       existence <- numerator > 0
                       logratios <- rep(0, length(numerator))
                       logratios[existence] <- log(numerator[existence] / refcounts[existence])
                       if (!all(existence)){
                         logratios[!existence] <- log(1/(refcounts[!existence]+1))
                       }
                       design_matrix <- cbind(1, metadata[, c(covar, adjust)])

                       LRT_result <- LRT_censored_regression(Y=logratios, Delta=existence, X=as.matrix(design_matrix),
                                                             dist=ratio_model, Firth=T)

                       beta_estimate <- LRT_result$full_estim_param[2] |> unname()
                       teststat_beta <- LRT_result$teststat
                       pval <- LRT_result$pval
                       estimation <- list(ID=j, effect=beta_estimate,
                                          teststat=teststat_beta, pval=pval)

                       return(as.data.frame(estimation))
                     }

  stopCluster(cl)
  CR_result$effect[outcome$ID] <- outcome$effect
  CR_result$teststat[outcome$ID] <- outcome$teststat
  CR_result$pval[outcome$ID] <- outcome$pval

  return(CR_result)
}
