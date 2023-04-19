
source('ANCOM+/overdisperse_GLM.R')
source('ANCOM+/Reference_Selection.R')

library(ACAT)
cauchy_combination <- function(pvalvec){
  pvalvec <- na.omit(pvalvec)
  return(ACAT(pvalvec))
}


polda <- function(otu_table, metadata, covar){

  taxa_names <- row.names(otu_table) # assume taxa are rows
  glm_result <- pairwise_GLM(otu_table, metadata, covar) # overdispersed GLM
  reftaxa <- backward_selection(otu_table, glm_result) # reference taxa selection

  # Set up a matrix of Logistic Regression P Values
  glm_pvals <- glm_result %>% select(T1, T2, pval)
  glm_pvals_wide <- reshape(glm_pvals, idvar="T1",
                            timevar="T2", direction="wide")

  cnames <- colnames(glm_pvals_wide)
  cnames <- sub("pval.", "", cnames)
  colnames(glm_pvals_wide) <- cnames

  pval_upper <- glm_pvals_wide[, -1] |> as.matrix()
  pval_mat <- matrix(NA, nrow=nrow(pval_upper)+1, ncol=ncol(pval_upper)+1)
  pval_mat[upper.tri(pval_mat)] <- pval_upper[upper.tri(pval_upper, diag=TRUE)]
  pval_mat[lower.tri(pval_mat)] <- 0
  pval_mat <- pval_mat + t(pval_mat)

  rownames(pval_mat) <- row.names(otu_table)
  colnames(pval_mat) <- row.names(otu_table)

  # only select a subset of logistic regression tests
  pval_mat_subset <- pval_mat[!taxa_names %in% reftaxa, reftaxa]
  Taxa_TBD <- row.names(pval_mat_subset)

  final_pvals <- apply(pval_mat_subset, 1, cauchy_combination)
  final_adjusted_pvals <- p.adjust(final_pvals, method="BH")
  pval_df <- data.frame(Taxa=Taxa_TBD, Pval=final_adjusted_pvals)

  DiffTaxa <- Taxa_TBD[final_adjusted_pvals < 0.05]

  result <- list(Reference_Taxa = reftaxa,
                 DA_taxa = DiffTaxa,
                 P_Value = pval_df)

  return(result)
}

