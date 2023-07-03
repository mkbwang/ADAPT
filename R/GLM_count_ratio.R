

#' Overdispersed and Bias-corrected Logistic Regression for Count Ratios
#'
#' @param count_data microbiome abundance matrix. All entries are integers
#' @param metadata sample metadata containing the covariate of interest and confounding factors for adjustment
#' @param covar the name of the covariate of interest
#' @param adjust the names of the variables for adjustment, NULL if no adjustments
#' @param reftaxa the list of reference taxa
#' @param complement If true, model the count ratios between each taxa beyond the reference set and the sum of taxa in the reference set. If false, model the relative abundance of taxa within the reference set.
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom stats pnorm
#' @returns The effect sizes, standard errors and p values for each taxon
#' @export
GLM_count_ratio <- function(count_data, metadata, covar, adjust=NULL, reftaxa=NULL, complement=FALSE){
  alltaxa <- rownames(count_data)
  if (is.null(reftaxa)) reftaxa <- alltaxa
  TBDtaxa <- NULL
  if (complement){ # examine count ratio between other taxa and the subset
    TBDtaxa <- setdiff(alltaxa, reftaxa)
  } else{ # relative abundance within the subset of taxa
    TBDtaxa <- reftaxa
  }
  refcounts <- colSums(count_data[reftaxa, ])
  TBD_counts <- count_data[TBDtaxa, ]
  num_taxa <- length(TBDtaxa)
  glmdisp_result <- data.frame(Taxon = TBDtaxa,
                               effect=0,
                               SE=0,
                               teststat=0,
                               pval=0)

  cores <- availableCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)

  outcome <- foreach(j=1:nrow(glmdisp_result), .combine=rbind,
                     .inorder=FALSE) %dopar% {

                       test_taxon <- as.character(glmdisp_result$Taxon[j])
                       numerator <- TBD_counts[test_taxon, ] # the counts of taxa of interest
                       denominator <-  numerator + refcounts
                       design_matrix <- cbind(1, metadata[, c(covar, adjust)])
                       # remove the samples with zero counts for both numerator and denominator
                       null_filter <- denominator != 0
                       numerator <- numerator[null_filter]
                       denominator <- denominator[null_filter]
                       design_matrix <- design_matrix[null_filter, ]

                       logit_result <- overdispersed_firth_IRWLS(y_vec=numerator, m_vec=denominator,
                                                                 X_mat=as.matrix(design_matrix))

                       beta_estimate <- logit_result$beta_hat[2]
                       se_beta_estimate <- sqrt(logit_result$vcov_beta[2, 2])
                       teststat_beta <- beta_estimate/se_beta_estimate
                       pval <- pnorm(-abs(teststat_beta))*2
                       estimation <- list(ID=j, effect=beta_estimate, SE=se_beta_estimate,
                                          teststat=teststat_beta, pval=pval)

                       return(as.data.frame(estimation))
                     }

  stopCluster(cl)
  glmdisp_result$effect[outcome$ID] <- outcome$effect
  glmdisp_result$SE[outcome$ID] <- outcome$SE
  glmdisp_result$teststat[outcome$ID] <- outcome$teststat
  glmdisp_result$pval[outcome$ID] <- outcome$pval

  return(glmdisp_result)
}
