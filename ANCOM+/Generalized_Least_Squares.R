## ----setup, include=FALSE----------------------------------------------------------
library(Matrix)
library(lme4)

# Woodbury inverse
woodbury_inverse <- function(diagA, U_mat, diagC, V_mat){
  # calculate (A + UCV)^-1, A and C being diagonal matrix
  A_inv <- sparseMatrix(i=1:length(diagA), j=1:length(diagA),
                        x=1/diagA) # feed in a vector, get a sparse matrix
  C_inv <- sparseMatrix(i=1:length(diagC), j=1:length(diagC),
                        x=1/diagC)

  C_VAU <- C_inv + V_mat %*% A_inv %*% U_mat
  C_VAU_inv <- chol2inv(chol(C_VAU))
  C_VAU_inv <- drop0(C_VAU_inv, tol=1e-12)
  output <- A_inv - A_inv %*% U_mat %*% C_VAU_inv %*% V_mat %*% A_inv
  output <- drop0(output, tol=1e-12)
  return(output)
}

# pseudo inverse
symmetric_pseudo_inverse <- function(symmat, rank){
  # Moore-Penrose inverse for positive semidefinite matrix with no full rank
  eigen_result <- eigen(symmat, symmetric = TRUE)
  eigen_values <- eigen_result$values
  eigen_vectors <- eigen_result$vectors
  inv_eig_values <- rep(0, length(eigen_values))
  inv_eig_values[1:rank] <- 1 / eigen_values[1:rank]
  inv_mat <- eigen_vectors %*% diag(inv_eig_values) %*% t(eigen_vectors)
  return(inv_mat)
}


lmer_function <- function(count_data, glm_result, reference=c("mean", "median")){

  num_taxa <- nrow(count_data)
  taxa_names <- rownames(count_data)

  # get the Z scores of all the logistic regressions
  Z_score <- glm_result$effect / glm_result$SE
  Z_score <- c(0, Z_score)

  # design matrix, contains one row for anchoring tau_1 at zero
  allpairs <- combn(seq(1, num_taxa), m=2)
  num_pair <- ncol(allpairs)
  design_matrix_rowid <- c(1, rep(seq(2, num_pair+1), 2))
  design_matrix_colid <- c(1, allpairs[1,], allpairs[2,])
  entry_vals <- c(1, rep(c(1, -1), each=num_pair))
  design_matrix <- sparseMatrix(i=design_matrix_rowid, j=design_matrix_colid,
                                x=entry_vals) |> as.matrix()

  data_df <- as.data.frame(design_matrix)
  colnames(data_df) <- taxa_names
  data_df$Z <- Z_score
  data_df$T1 <- c(glm_result$T1[1], glm_result$T1) |> as.factor()

  formula_str <- sprintf("Z ~ %s - 1 + (1|T1) ",
                         paste(taxa_names, collapse="+"))
  regression_formula <- as.formula(formula_str)
  model_fitting <- lmer(regression_formula, data=data_df)
  model_summary <- summary(model_fitting)
  tau_hat <- model_summary$coefficients[, 1]
  median_tau_hat <- median(tau_hat)
  mean_tau_hat <- mean(tau_hat)
  if (reference == "mean"){
  shifted_tau_hat <- tau_hat - mean_tau_hat
  } else{
   shifted_tau_hat <- tau_hat - median_tau_hat
  }
  stderror_tau <- model_summary$coefficients[, 2]
  tau_test_statistic <- shifted_tau_hat / stderror_tau
  pval <- pnorm(-abs(tau_test_statistic)) * 2
  pval_adjusted <- p.adjust(pval, method="BH")

  result <- data.frame(Taxa = taxa_names,
                      tau_hat = shifted_tau_hat,
                      stderror_tau = stderror_tau,
                      tau_test_statistic=tau_test_statistic,
                      Pval=pval,
                      DA=pval < 0.05)
  return(list(result, model_summary))
}



# generalized least square
generalized_ls <- function(count_data, metadata, glm_result, reference=c("mean", "median"), corr=TRUE){
  num_taxa <- nrow(count_data)

  # get the Z scores of all the logistic regressions
  Z_score <- glm_result$effect / glm_result$SE

  # design matrix
  allpairs <- combn(seq(1, num_taxa), m=2)
  num_pair <- ncol(allpairs)
  design_matrix_rowid <- rep(seq(1, num_pair), 2)
  design_matrix_colid <- c(allpairs[1,], allpairs[2,])
  entry_vals <- rep(c(1, -1), each=num_pair)
  design_matrix <- sparseMatrix(i=design_matrix_rowid, j=design_matrix_colid,
                                x=entry_vals)

  dm_rank <- rankMatrix(design_matrix, method='qr')[[1]]

  nonzero_proportion <- rowMeans(count_data != 0)
  taxa_reffect_var <- rep(1, num_taxa) # random effect variance
  residual_var <- rep(1, num_pair) # residual variance

  # covariance matrix inverse(without the unknown sigma^2 factor)
  V_inv <- woodbury_inverse(residual_var, design_matrix, taxa_reffect_var,
                            t(design_matrix))

  XtVX <- t(design_matrix) %*% V_inv %*% design_matrix
  XtVX <- as.matrix(XtVX)
  XtVX_inv <- symmetric_pseudo_inverse(XtVX, rank=dm_rank)

  XtX <- t(design_matrix) %*% design_matrix
  XtX <- as.matrix(XtX)
  XtX_inv <- symmetric_pseudo_inverse(XtX, rank=dm_rank)

  if (corr){
    # estimate tau
    tau_hat <- XtVX_inv %*% t(design_matrix) %*% V_inv %*% Z_score
    tau_hat <- as.vector(tau_hat)
    # residual
    resid <- Z_score - design_matrix %*% tau_hat
    # sigma^2 estimate (MLE)
    sigma2_hat <- t(resid) %*% V_inv %*% resid / nrow(resid)
    # estimated covariance of tau
    covar_tau <- drop(sigma2_hat) * XtVX_inv
    stderror_tau <- sqrt(diag(covar_tau))
  } else{
    # estimate tau
    tau_hat <- XtX_inv %*% t(design_matrix) %*% Z_score
    tau_hat <- as.vector(tau_hat)
    # residual
    resid <- Z_score - design_matrix %*% tau_hat
    # sigma^2 estimate (MLE)
    sigma2_hat <- t(resid)  %*% resid / nrow(resid)
    # estimated covariance of tau
    covar_tau <- drop(sigma2_hat) * XtX_inv
    stderror_tau <- sqrt(diag(covar_tau))
  }

  if (reference == "median"){
    shift <- median(tau_hat) # shift the estimated taus so that the median is zero
    shifted_tau_hat <- tau_hat - shift
  } else{ # if mean, no need to shift
    shifted_tau_hat <- tau_hat
  }

  test_statistic <- as.vector(shifted_tau_hat / stderror_tau)
  pval <- 2 * pnorm(-abs(test_statistic))

  adjusted_pval <- p.adjust(pval, method="BH")
  result <- data.frame(Taxa_name = rownames(count_data),
                       tau_hat = shifted_tau_hat,
                       sderror = stderror_tau,
                       pval_raw = pval,
                       pval_adjust = adjusted_pval,
                       DA = adjusted_pval < 0.05)

  return(result)
}

