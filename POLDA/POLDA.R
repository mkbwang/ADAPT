
source('POLDA/overdisperse_GLM.R')
source('POLDA/Reference_Selection.R')

library(ACAT)
cauchy_combination <- function(pvalvec){
  pvalvec <- na.omit(pvalvec)
  return(ACAT(pvalvec))
}


polda <- function(otu_table, metadata, covar, covartype=c("categorical", "numerical")){

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

  glm_result <- pairwise_GLM(otu_table, metadata, covar) # overdispersed GLM
  reftaxa <- backward_selection(otu_table, glm_result) # reference taxa selection

  # Set up a matrix of Logistic Regression P Values
  glm_pvals <- glm_result %>% dplyr::select(T1, T2, pval)
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

  result <- list(Structural_zero_Taxa = struct_zero_taxa,
                 Reference_Taxa = reftaxa,
                 DA_taxa = DiffTaxa,
                 P_Value = pval_df)

  return(result)
}

