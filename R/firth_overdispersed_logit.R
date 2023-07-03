
sigmoid <- function(eta_vec){
  exp(eta_vec)/(1+exp(eta_vec))
}



vanilla_IRWLS <- function(y_vec, m_vec, X_mat, phi=0, tol=1e-4){
  beta_vec <- rep(0, ncol(X_mat))
  weights <- 1/(1 + (m_vec - 1) * phi)
  XWX <- NULL
  XWX_inv <- NULL
  while(1){
    eta_vec <- as.vector(X_mat %*% beta_vec) # link
    pi_vec <- sigmoid(eta_vec) # estimated probability
    W <- diag(m_vec * pi_vec * (1-pi_vec) * weights)
    XWX <- t(X_mat) %*% W %*% X_mat # information matrix
    XWz <- XWX %*% beta_vec + t(X_mat) %*% (y_vec*weights - m_vec*pi_vec*weights)
    new_beta_vec <- solve(XWX, XWz)
    if (max(abs(beta_vec - new_beta_vec)) < tol) break
    beta_vec <- new_beta_vec
  }
  XWX_inv <- chol2inv(chol(XWX))
  return(list(beta_hat=beta_vec, vcov_beta=XWX_inv))
}


firth_IRWLS <- function(y_vec, m_vec, X_mat, phi=0, tol=1e-4){
  beta_vec <- rep(0, ncol(X_mat))
  weights <- 1/(1 + (m_vec - 1) * phi)
  XWX <- NULL
  XWX_inv <- NULL
  while(1){
    eta_vec <- as.vector(X_mat %*% beta_vec) # link
    pi_vec <- sigmoid(eta_vec) # estimated probability
    W <- diag(m_vec * pi_vec * (1-pi_vec) * weights)
    XWX <- t(X_mat) %*% W %*% X_mat # information matrix
    XWX_inv <- chol2inv(chol(XWX))
    H_mat <- X_mat %*% XWX_inv %*% t(X_mat) %*% W
    h_diag <- diag(H_mat)
    XWz <- XWX %*% beta_vec + t(X_mat) %*% (y_vec*weights - m_vec*pi_vec*weights)+
      t(X_mat) %*% (h_diag * (0.5 - pi_vec))
    new_beta_vec <- solve(XWX, XWz)
    if (max(abs(beta_vec - new_beta_vec)) < tol) break
    beta_vec <- new_beta_vec
  }
  return(list(beta_hat=beta_vec, vcov_beta=XWX_inv))
}


overdispersed_vanilla_IRWLS <- function(y_vec, m_vec, X_mat, epsilon=1e-5){
  phi <- 0  # initially no dispersion
  LRresult <- vanilla_IRWLS(y_vec=y_vec, m_vec=m_vec, X_mat=X_mat, phi=phi) # fit standard LR
  nobs <- nrow(X_mat)
  nfeature <- ncol(X_mat)
  beta_vec <- LRresult$beta_hat
  weights <- 1/((m_vec - 1) * phi + 1)
  XWX <- NULL
  XWX_inv <- NULL
  while(1){
    eta_vec <- as.vector(X_mat %*% beta_vec) # link
    pi_vec <- sigmoid(eta_vec) # estimated probability
    chi_2 <- sum((y_vec - m_vec*pi_vec)^2 * weights / (m_vec * pi_vec * (1-pi_vec))) # pearson statistic
    W <- diag(m_vec * pi_vec * (1-pi_vec) * weights)
    XWX <- t(X_mat) %*% W %*% X_mat # information matrix
    XWX_inv <- chol2inv(chol(XWX))
    if (chi_2 / (nobs-nfeature) - 1 < epsilon) break # no dispersion any more
    H_mat <- X_mat %*% XWX_inv %*% t(X_mat) %*% W # hat matrix
    hdiag <- diag(H_mat)
    # update overdispersion parameter
    phi <- (chi_2 - sum(weights * (1-hdiag))) / sum(weights * (m_vec - 1) * (1 - hdiag))
    if (phi < 0) break
    ## refit the model with the new phi
    weights <- 1/((m_vec - 1) * phi + 1)
    LRresult <- vanilla_IRWLS(y_vec, m_vec, X_mat, phi=phi, tol=1e-4)
    beta_vec <- LRresult$beta_hat
  }
  return(list(beta_hat=beta_vec, vcov_beta=XWX_inv, phi=phi))
}



overdispersed_firth_IRWLS <- function(y_vec, m_vec, X_mat, epsilon=1e-2){
  phi <- 0  # initially no dispersion
  LRresult <- firth_IRWLS(y_vec=y_vec, m_vec=m_vec, X_mat=X_mat, phi=phi) # fit standard LR
  nobs <- nrow(X_mat)
  nfeature <- ncol(X_mat)
  beta_vec <- LRresult$beta_hat
  weights <- 1/((m_vec - 1) * phi + 1)
  XWX <- NULL
  XWX_inv <- NULL
  while(1){
    eta_vec <- as.vector(X_mat %*% beta_vec) # link
    pi_vec <- sigmoid(eta_vec) # estimated probability
    chi_2 <- sum((y_vec - m_vec*pi_vec)^2 * weights / (m_vec * pi_vec * (1-pi_vec))) # pearson statistic
    W <- diag(m_vec * pi_vec * (1-pi_vec) * weights)
    XWX <- t(X_mat) %*% W %*% X_mat # information matrix
    XWX_inv <- chol2inv(chol(XWX))
    if (chi_2 / (nobs-nfeature) - 1 < epsilon) break # no dispersion any more
    H_mat <- X_mat %*% XWX_inv %*% t(X_mat) %*% W # hat matrix
    hdiag <- diag(H_mat)
    # update overdispersion parameter
    phi <- (chi_2 - sum(weights * (1-hdiag))) / sum(weights * (m_vec - 1) * (1 - hdiag))
    if (phi < 0) break
    ## refit the model with the new phi
    weights <- 1/((m_vec - 1) * phi + 1)
    LRresult <- firth_IRWLS(y_vec, m_vec, X_mat, phi=phi)
    beta_vec <- LRresult$beta_hat
  }
  return(list(beta_hat=beta_vec, vcov_beta=XWX_inv, phi=phi))
}

