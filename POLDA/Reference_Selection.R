## ----setup, include=FALSE----------------------------------------------------------
library(Matrix)

# pseudo inverse
symmetric_pseudo_inverse <- function(symmat, matrank){
  # Moore-Penrose inverse for positive semidefinite matrix with no full rank
  eigen_result <- eigen(symmat, symmetric = TRUE)
  eigen_values <- eigen_result$values
  eigen_vectors <- eigen_result$vectors
  inv_eig_values <- rep(0, length(eigen_values))
  inv_eig_values[1:matrank] <- 1 / eigen_values[1:matrank]
  inv_mat <- eigen_vectors %*% diag(inv_eig_values) %*% t(eigen_vectors)
  return(inv_mat)
}


# weighted least squares
weighted_ls <- function(design_mat, weights, Y, matrank, zero_counts, start=c("median", "mean", "mode")){

  start <- match.arg(start) # the first taxon to drop: median or mean?
  # weighted least squares
  XtWX <- t(design_mat) %*% weights %*% design_mat
  XtWX <- as.matrix(XtWX)
  if (matrank < ncol(design_mat)){
    XtWX_inv <- symmetric_pseudo_inverse(XtWX, matrank=matrank)
  } else{
    XtWX_inv <- chol2inv(chol(XtWX))
  }
  # estimate tau
  tau_hat <- XtWX_inv %*% t(design_mat) %*% weights %*% Y
  tau_hat <- as.vector(tau_hat)
  # residual
  resid <- Y - design_mat %*% tau_hat
  # sigma^2 estimate (MLE)
  sigma2_hat <- t(resid) %*% weights %*% resid / nrow(resid)
  # estimated covariance of tau
  covar_tau <- drop(sigma2_hat) * XtWX_inv
  stderror_tau <- sqrt(diag(covar_tau))
  # find the next variable to drop
  test_statistic <- as.vector(tau_hat / stderror_tau)
  if (matrank < ncol(design_mat) & start == "median"){
    # drop the median at first
    rank2median <- rank(abs(tau_hat-median(tau_hat)))
    revised_rank <- rank2median + zero_counts
    drop_id <- which.min(revised_rank)
  } else if (matrank < ncol(design_mat) & start == "mode"){
    kernel_fit_tau <- density(tau_hat, bw="SJ", kernel="epanechnikov")
    tau_mode <- kernel_fit_tau$x[which.max(kernel_fit_tau$y)]
    rank2mode <- rank(abs(tau_hat - tau_mode))
    revised_rank <- rank2mode + zero_counts
    drop_id <- which.min(revised_rank)
  } else{
    original_rank <- rank(abs(test_statistic))
    revised_rank <- original_rank + zero_counts
    drop_id <- which.min(revised_rank)
  }
  # AIC calculation
  weighted_RSS <- drop(t(resid) %*% weights %*% resid)
  AIC_value <- nrow(design_mat) * log(weighted_RSS) +
    ncol(design_mat) * 2
  return(list(AIC=AIC_value, drop=drop_id, tau_hat=tau_hat))
}


# generalized least square
backward_selection <- function(count_data, glm_result, start=c("median", "mean", "mode")){

  start <- match.arg(start) # the first taxon to drop: median or mean?
  num_taxa <- nrow(count_data)
  taxa_names <- rownames(count_data)
  zero_counts <- rowSums(count_data == 0)
  # some logistic regression may fail due to quasi separation
  success_pairs <- !is.na(glm_result$SE)
  estimated_effect <- glm_result$effect
  estimated_variance <- (glm_result$SE)^2

  # design matrix
  allpairs <- combn(seq(1, num_taxa), m=2)

  num_pair <- ncol(allpairs)
  design_matrix_rowid <- rep(seq(1, num_pair), 2)
  design_matrix_colid <- c(allpairs[1, ], allpairs[2, ])
  entry_vals <- rep(c(1, -1), each=num_pair)
  design_matrix <- sparseMatrix(i=design_matrix_rowid, j=design_matrix_colid,
                                x=entry_vals)

  design_matrix <- design_matrix[success_pairs, ]
  # covariance matrix inverse(without the unknown sigma^2 factor) as weights
  weights <- sparseMatrix(i=seq(1, sum(success_pairs)), j=seq(1, sum(success_pairs)),
                          x=1 / estimated_variance[success_pairs])
  logit_effects <- estimated_effect[success_pairs]

  # first fit the full model
  current_model <- weighted_ls(design_matrix, weights, logit_effects, num_taxa-1, 
                               zero_counts=zero_counts, start)
  full_tau_hat <- current_model$tau_hat
  full_AIC <- current_model$AIC
  names(full_tau_hat) <- taxa_names
  AIC_path <- full_AIC # AIC trajectory
  reference_taxa <- c()
  # backward selection
  while(TRUE){
    current_AIC <- current_model$AIC
    dropid <- current_model$drop
    design_matrix <- design_matrix[, -dropid]
    zero_counts <- zero_counts[-dropid]
    new_model <- weighted_ls(design_matrix, weights, logit_effects, 
                             zero_counts=zero_counts, ncol(design_matrix))
    new_AIC <- new_model$AIC
    AIC_path <- c(AIC_path, new_AIC)
    if (new_AIC < full_AIC){
      reference_taxa <- c(reference_taxa, taxa_names[dropid])
      taxa_names <- taxa_names[-dropid]
      current_model <- new_model
    } else{
      break
    }
  }

  num_reftaxa <- which.min(AIC_path) - 1
  return(list(reftaxa = reference_taxa[1:num_reftaxa],
              AIC_traj = AIC_path,
              full_tau_hat = full_tau_hat))
}

