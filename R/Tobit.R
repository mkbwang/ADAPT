


## hessian(negative information) matrix of Gaussian distribution
gaussian_hessian <- function(theta, Y, Delta, X, fixed=NULL){
  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  tau <- theta[1:(num_params-1)] # coefficients for location
  if(!is.null(fixed)){
    tau[fixed] <- 0
  }
  phi <- exp(theta[num_params]) # inverse scale parameter

  z <- phi*Y - as.vector(X %*% tau)

  cumz <- pnorm(z)# cumulative gaussian distribution of z
  exp_z2 <- exp(-z^2)

  deriv1_z <- -Delta*z + (1-Delta) / cumz * sqrt(exp_z2) / sqrt(2*pi)
  deriv2_z <- -Delta - (1-Delta)/(2*pi)*exp_z2/(cumz^2) - (1-Delta)/sqrt(2*pi)*z*sqrt(exp_z2)/cumz

  output_hessian <- matrix(0, nrow=num_params, ncol=num_params)
  ## first calculate the hessian for tau only
  output_hessian[1:(num_params-1), 1:(num_params-1)] <- t(X) %*% diag(deriv2_z) %*% X

  ## then calculate the hessian for phi
  output_hessian[num_params, num_params] <- sum(deriv2_z * Y^2) + sum(-Delta/(phi^2))

  ## finally calculate the derivative over both tau and phi
  weights <- - deriv2_z * Y
  output_hessian[num_params, 1:(num_params-1)] <- output_hessian[1:(num_params-1), num_params] <-
    as.vector(weights %*% X)

  return(output_hessian)
}

## gaussian distribution log likelihood
gaussian_llk <- function(theta, Y, Delta, X, Firth=T, fixed=NULL, pen=0.5){
  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  tau <- theta[1:(num_params-1)] # coefficients for location
  if (!is.null(fixed)){
    tau[fixed] <- 0
  }
  phi <- exp(theta[num_params]) # scale parameter
  z <- phi*Y - as.vector(X %*% tau)

  llk <- sum(Delta * (-0.5 * log(2 * pi) + log(phi) - 0.5 * z^2)) + sum((1-Delta) * log(pnorm(z)))

  if(Firth){ # add penalty
    info_mat <- -gaussian_hessian(theta, Y, Delta, X, fixed=fixed)
    llk <- llk + pen*(determinant(info_mat, logarithm=T)$modulus[1])
  }

  return(llk)
}



## hessian(negative information) matrix of logistic distribution
logistic_hessian <- function(theta, Y, Delta, X, fixed=NULL){
  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  tau <- theta[1:(num_params-1)] # coefficients for location
  if(!is.null(fixed)){
    tau[fixed] <- 0
  }
  phi <- exp(theta[num_params]) # inverse scale parameter

  z <- phi*Y - as.vector(X %*% tau)

  deriv1_z <- -Delta + (1+Delta) / (1+exp(z))
  deriv2_z <- - (1+Delta) * exp(z)/((1+exp(z))^2)

  output_hessian <- matrix(0, nrow=num_params, ncol=num_params)
  ## first calculate the hessian for tau only
  output_hessian[1:(num_params-1), 1:(num_params-1)] <- t(X) %*% diag(deriv2_z) %*% X

  ## then calculate the hessian for phi
  output_hessian[num_params, num_params] <- sum(deriv2_z * Y^2) + sum(-Delta/(phi^2))

  ## finally calculate the derivative over both tau and phi
  weights <- - deriv2_z * Y
  output_hessian[num_params, 1:(num_params-1)] <- output_hessian[1:(num_params-1), num_params] <-
    as.vector(weights %*% X)

  return(output_hessian)
}


## logistic distribution log likelihood
logistic_llk <- function(theta, Y, Delta, X, Firth=T, fixed=NULL, pen=0.5){
  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  tau <- theta[1:(num_params-1)] # coefficients for location
  if (!is.null(fixed)){
    tau[fixed] <- 0
  }
  phi <- exp(theta[num_params]) # scale parameter
  z <- phi*Y - as.vector(X %*% tau)

  llk <- sum(Delta*(log(phi)-z)) - sum((1+Delta)*log(1+exp(-z)))

  if(Firth){ # add penalty
    info_mat <- -logistic_hessian(theta, Y, Delta, X, fixed=fixed)
    llk <- llk + pen*(determinant(info_mat, logarithm=T)$modulus[1])
  }

  return(llk)
}



# estimate censored regression coefficients
estim_censored_regression <- function(Y, Delta, X, Firth=T, fixed=NULL, dist=c("lognormal", "loglogistic"), pen=0.5){

  dist <- match.arg(dist)

  selected_llk <- NULL # likelihood function
  initial_intercept <- mean(Y)/sqrt(var(Y))
  initial_tau <- rep(0, ncol(X)-1)
  initial_logphi <- log(1/sqrt(var(Y)))

  selected_llk <- NULL
  if (dist == "lognormal"){
    selected_llk <- gaussian_llk
  } else{ # loglogistic
    selected_llk <- logistic_llk
  }

  initial_param <- c(initial_intercept, initial_tau, initial_logphi)
  num_params <- length(initial_param)
  lower_bound <- c(rep(-Inf, num_params-1), 0)
  upper_bound <- rep(Inf, num_params)
  estim_params <- NULL

  estim_result <- optim(par=initial_param, fn=selected_llk,
                        Y=Y, Delta=Delta, X=X, Firth=Firth, fixed=fixed, pen=pen,
                        control=list(fnscale=-1), method="BFGS")
  estim_params <- estim_result$par


  return(estim_params)
}


LRT_censored_regression <- function(Y, Delta, X,  Firth=T, dist=c("lognormal", "loglogistic"), test_param=2, pen=0.5){


  dist <- match.arg(dist)

  full_estim_param <- NULL
  reduced_estim_param <- NULL
  teststat <- NULL
  pval <- NULL

  selected_llk <- NULL
  if (dist == "lognormal"){
    selected_llk <- gaussian_llk
  } else{ # loglogistic
    selected_llk <- logistic_llk
  }


  full_estim_param <- estim_censored_regression(Y=Y, Delta=Delta, X=X, Firth=Firth, dist=dist,
                                                fixed=NULL, pen=pen)
  full_llk <- selected_llk(theta=full_estim_param,
                           Y=Y, Delta=Delta, X=X, Firth=Firth, pen=pen)
  reduced_estim_param <- estim_censored_regression(Y=Y, Delta=Delta, X=X, Firth=Firth, dist=dist,
                                                   fixed=test_param, pen=pen)
  reduced_llk <- selected_llk(theta=reduced_estim_param,
                              Y=Y, Delta=Delta, X=X, Firth=Firth, pen=pen)

  teststat <- 2*(full_llk - reduced_llk)
  pval <- 1-pchisq(teststat, df=length(test_param))

  num_param <- length(full_estim_param)
  # scale parameter
  full_estim_param[num_param] <- 1/exp(full_estim_param[num_param])
  reduced_estim_param[num_param] <- 1/exp(reduced_estim_param[num_param])
  # effect sizes
  full_estim_param[1:(num_param-1)] <- full_estim_param[1:(num_param-1)] * full_estim_param[num_param]
  reduced_estim_param[1:(num_param-1)] <- reduced_estim_param[1:(num_param-1)] * reduced_estim_param[num_param]

  names(full_estim_param)[1:(num_param-1)] <- colnames(X)
  names(full_estim_param)[num_param] <- "Scale"
  names(reduced_estim_param)[1:(num_param-1)] <- colnames(X)
  names(reduced_estim_param)[num_param] <- "Scale"
  result <- list(full_estim_param=full_estim_param,
                 reduced_estim_param=reduced_estim_param,
                 teststat=teststat,
                 pval=pval)

  return(result)
}


# library(survival)
#
# #test data
# set.seed(2035)
# time1 <- rexp(n=20, rate=0.1)
# time2 <- rexp(n=20, rate=0.1)
# alltimes <- c(time1, time2)
# event <- sample(c(0, 1), size=40, replace=TRUE, prob=c(0.8, 0.2))
# alltimes[!event] <- alltimes[!event] / 2
# neglogtime <- -log(alltimes)
# covariate <- c(rep(0, 20), rep(1, 20))
# adjust_covar <- rnorm(n=40)
# Xmat <- cbind(1, covariate)
#
#
#
# # AFT model with log normal regression
# lognormal_vanilla_result <- LRT_censored_regression(Y=neglogtime, Delta=event, X=Xmat, dist="lognormal",
#                                                     Firth=F)
#
#
# standard_full <- survreg(Surv(neglogtime, event, type="left") ~ covariate, dist="gaussian")
# standard_null <- survreg(Surv(neglogtime, event, type="left") ~ 1, dist="gaussian")
#
# lognormal_firth_result <- LRT_censored_regression(Y=neglogtime, Delta=event, X=Xmat, dist="lognormal",
#                                                   Firth=T, pen=0.5)
#
#
# loglogistic_vanilla_result <- LRT_censored_regression(Y=neglogtime, Delta=event, X=Xmat, dist="loglogistic",
#                                                       Firth=F)
#
#
# standard_full <- survreg(Surv(neglogtime, event, type="left") ~ covariate, dist="logistic")
# standard_null <- survreg(Surv(neglogtime, event, type="left") ~ 1, dist="logistic")
#
#
# loglogistic_firth_result <- LRT_censored_regression(Y=neglogtime, Delta=event, X=Xmat, dist="loglogistic",
#                                                       Firth=T, pen=0.5)


# info_mat <- -gaussian_hessian(lognormal_vanilla_result$full_estim_param, Y=neglogtime,
#                               Delta=event, X=Xmat, fixed=NULL)
#
# full_llk <- gaussian_llk(theta=lognormal_vanilla_result$full_estim_param,
#                          Y=neglogtime,
#                          Delta=event, X=Xmat, Firth=F)
#
# full_llk_pen <- gaussian_llk(theta=lognormal_firth_result$full_estim_param,
#                              Y=neglogtime,
#                              Delta=event, X=Xmat, Firth=T)
# reduced_llk <- gaussian_llk(theta=lognormal_vanilla_result$reduced_estim_param,
#                             Y=neglogtime,
#                             Delta=event, X=Xmat, Firth=F)
#
# reduced_llk_pen <- gaussian_llk(theta=lognormal_firth_result$reduced_estim_param,
#                                 Y=neglogtime,
#                                 Delta=event, X=Xmat, Firth=T)
#
# info_mat_reduced <- -gaussian_hessian(lognormal_firth_result$reduced_estim_param, Y=neglogtime,
#                                       Delta=event, X=Xmat)
#
# determinant(info_mat, logarithm=T)$modulus[1]
# determinant(info_mat_reduced, logarithm=T)$modulus[1]


