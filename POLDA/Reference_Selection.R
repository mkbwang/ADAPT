## ----setup, include=FALSE----------------------------------------------------------
library(Matrix)

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


# weighted least squares
weighted_ls <- function(design_mat, weights, Y, dropids=c()){
  
  if (length(dropids) > 0){
    design_mat <- design_mat[, -dropids]    
  }
  # weighted least squares
  XtWX <- t(design_mat) %*% weights %*% design_mat
  XtWX <- as.matrix(XtWX)
  if (length(dropids) == 0){ # no drop, design matrix not full rank
    XtWX_inv <- symmetric_pseudo_inverse(XtWX, rank=ncol(design_mat) - 1)
  } else{
    XtWX_inv <- chol2inv(chol(XtWX))
  }
  # estimate tau
  tau_hat <- XtWX_inv %*% t(design_mat) %*% weights %*% Y
  tau_hat <- as.vector(tau_hat)
  # residual
  resid <- Y - design_mat %*% tau_hat
  # # sigma^2 estimate (MLE)
  # sigma2_hat <- t(resid) %*% weights %*% resid / nrow(resid)
  # # estimated covariance of tau
  # covar_tau <- drop(sigma2_hat) * XtWX_inv
  # stderror_tau <- sqrt(diag(covar_tau))
  # # find the next variable to drop
  # test_statistic <- as.vector(tau_hat / stderror_tau)
  # if (rank < ncol(design_mat) & start == "median"){
  #   # drop the median at first
  #   drop_id <- which.min(abs(tau_hat-median(tau_hat)))
  # } else if (rank < ncol(design_mat) & start == "mode"){
  #   kernel_fit_tau <- density(tau_hat, bw="SJ", kernel="epanechnikov")
  #   tau_mode <- kernel_fit_tau$x[which.max(kernel_fit_tau$y)]
  #   drop_id <- which.min(abs(tau_hat - tau_mode))
  # } else{
  #   drop_id <- which.min(abs(test_statistic))
  # }
  # BIC calculation
  weighted_RSS <- drop(t(resid) %*% weights %*% resid)
  BIC_value <- nrow(design_mat) * log(weighted_RSS) +
    ncol(design_mat) * log(nrow(design_mat))
  return(list(BIC=BIC_value, tau_hat=tau_hat, dropids=dropids))
}


# generalized least square
backward_selection <- function(count_data, glm_result, start=c("median", "mean")){

  start <- match.arg(start) # the first taxon to drop: median or mean?
  num_taxa <- nrow(count_data)
  taxa_names <- rownames(count_data)
  
  # retrieve prevalence and abundance information and generate ranking based on basic stats
  zero_count <- rowSums(count_data == 0)
  # sequencing_depths <- colSums(count_data)
  # relative_abundance <- count_data / 
  #   t(replicate(nrow(count_data), sequencing_depths))
  # avg_rel_abundance <- rowMeans(relative_abundance) # average of relative abundance
  # var_rel_abundance <- apply(relative_abundance, 1, var) # variance of relative abundance
  # taxa_basic_stats <- cbind(prevalence, avg_rel_abundance, var_rel_abundance) |> as.data.frame()
  # taxa_basic_stats$Taxa_name <- row.names(taxa_basic_stats)
  # taxa_basic_stats  <- taxa_basic_stats %>% arrange(desc(prevalence),
  #                                                   var_rel_abundance,
  #                                                   desc(avg_rel_abundance))
  # taxa_basic_stats$basic_ranking <- seq(1, nrow(taxa_basic_stats))
  
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
  current_model <- weighted_ls(design_matrix, weights, logit_effects, dropids = c())
  full_tau_hat <- current_model$tau_hat
  names(full_tau_hat) <- taxa_names
  
  # rank taxa based on distance to mean/median of tau
  tau_ranking <- rank(abs(full_tau_hat)) # based on mean
  if (start == "median"){
    median_tau <- median(full_tau_hat)
    tau_ranking <- rank(abs(full_tau_hat - median_tau))
  }
  tau_ranking_df <- data.frame(Taxa_name = taxa_names,
                               Tau = full_tau_hat,
                               Tau_rank = tau_ranking)
  tau_ranking_df$Zero_Count <- zero_count
  
  # combine the ranking of tau and basic statistics
  tau_ranking_df$Tau_rank_adjusted <- tau_ranking_df$Tau_rank + tau_ranking_df$Zero_Count
  tau_ranking_df <- tau_ranking_df %>% arrange(Tau_rank_adjusted)
  # order of design matrix column to drop
  drop_order <- match(tau_ranking_df$Taxa_name, taxa_names)
  full_BIC <- current_model$BIC
  num_reftaxa <- 0
  BIC_values <- c()
  BIC_values[1] <- full_BIC
  # backward selection
  while(TRUE){
    current_BIC <- current_model$BIC
    new_model <- weighted_ls(design_matrix, weights, logit_effects, 
                             dropids=drop_order[seq(1, num_reftaxa+1)])
    new_BIC <- new_model$BIC
    BIC_values <- c(BIC_values, new_BIC)
    if (new_BIC < full_BIC){
      num_reftaxa <- num_reftaxa + 1
      current_model <- new_model
    } else{
      break
    }
  }

  return(list(reftaxa = tau_ranking_df$Taxa_name[seq(1, num_reftaxa)],
              taxa_info = tau_ranking_df,
              BIC_values = BIC_values))
}

