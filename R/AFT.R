
expit <- function(x){

  output <- rep(0, length(x))
  output[x >= 0] <- 1/(1+exp(-x[x >= 0]))
  output[x < 0] <- exp(x[x < 0])/(1+exp(x[x < 0]))

  return(output)
}

logoneplusexp <- function(x){

  output <- rep(0, length(x))
  smallabs <- abs(x) <= 30
  verylarge <- x > 30
  verysmall <- x < -30
  output[smallabs] <- log(1+exp(x[smallabs]))
  output[verylarge] <- x[verylarge]
  output[verysmall] <- 0

  return(output)
}



## hessian(negative information) matrix of logistic distribution
logistic_hessian <- function(theta, Y, Delta, X, fixed=NULL){
  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  if(!is.null(fixed)) {
    beta[fixed] <- 0
  }
  lambda <- theta[num_params] # scale parameter
  eta <- as.vector(X %*% beta)
  z <- (eta - Y)/exp(lambda)

  expit_z <- expit(z)

  output_hessian <- matrix(0, nrow=num_params, ncol=num_params)

  ## first calculate the hessian for the beta only
  weights <- -(1+Delta)/exp(2*lambda)*expit_z*(1-expit_z)
  output_hessian[1:(num_params-1), 1:(num_params-1)] <- t(X) %*% diag(weights) %*% X

  ## then calculate the hessian for lambda
  output_hessian[num_params, num_params] <- -sum(((1+Delta)*expit_z - Delta)*z) -
    sum((1+Delta)*expit_z*(1-expit_z)*z^2)

  ## finally calculate derivative over both beta and lambda
  weights <- 1/exp(lambda)*(1-(1+Delta)*(1-expit_z) + (1+Delta)*z*expit_z*(1-expit_z))
  output_hessian[1:(num_params-1), num_params] <- output_hessian[num_params, 1:(num_params-1)] <-
    as.vector(weights %*% X)

  return(output_hessian)
}

## logistic distribution log likelihood
logistic_llk <- function(theta, Y, Delta, X, Firth=T, fixed=NULL){

  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  if(!is.null(fixed)) {
    beta[fixed] <- 0
  }
  lambda <- theta[num_params] # scale parameter
  eta <- as.vector(X %*% beta)
  z <- (eta - Y)/exp(lambda)
  logonepluxexp_z <- sapply(z, logoneplusexp)
  llk <- sum(Delta * z)- sum(Delta * lambda) -
    sum((1+Delta) * logoneplusexp(z))

  if(Firth){ # add penalty
    info_mat <- -logistic_hessian(theta, Y, Delta, X, fixed=fixed)
    llk <- llk + 0.5*(determinant(info_mat, logarithm=T)$modulus[1])
  }

  return(llk)
}


## hessian(negative information) matrix of gumbel distribution
gumbel_hessian <- function(theta, Y, Delta, X, fixed=NULL){

  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  if(!is.null(fixed)){
    beta[fixed] <- 0
  }
  lambda <- theta[num_params] # scale parameter
  eta <- as.vector(X %*% beta)
  z <- (eta - Y)/exp(lambda)

  output_hessian <- matrix(0, nrow=num_params, ncol=num_params)
  ## first calculate the hessian for the beta only
  weights <- -exp(z-2*lambda)
  output_hessian[1:(num_params-1), 1:(num_params-1)] <- t(X) %*% diag(weights) %*% X


  ## then calculate the hessian for lambda only
  output_hessian[num_params, num_params] <- sum(Delta * z - exp(z)*z - exp(z)*z^2)

  ## finally calculate the derivative over both beta and lambda
  weights <- exp(z-lambda)*z+exp(z-lambda)-Delta*exp(-lambda)
  output_hessian[num_params, 1:(num_params-1)] <- output_hessian[1:(num_params-1), num_params] <-
    as.vector(weights %*% X)

  return(output_hessian)

}

## gumbel distribution log likelihood
gumbel_llk <- function(theta, Y, Delta, X, Firth=T, fixed=NULL){

  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  if (!is.null(fixed)){
    beta[fixed] <- 0
  }
  lambda <- theta[num_params] # scale parameter
  eta <- as.vector(X %*% beta)
  z <- (eta - Y)/exp(lambda)

  llk <- -sum(Delta * lambda) + sum(Delta * z) - sum(exp(z))

  if(Firth){ # add penalty
    info_mat <- -gumbel_hessian(theta, Y, Delta, X, fixed=fixed)
    llk <- llk + 0.5*(determinant(info_mat, logarithm=T)$modulus[1])
  }

  return(llk)

}

# estimate censored regression coefficients
estim_censored_regression <- function(Y, Delta, X, Firth=T, fixed=NULL, dist=c("loglogistic", "weibull")){
  selected_llk <- NULL # likelihood function
  initial_intercept <- mean(Y)
  initial_beta <- rep(0, ncol(X)-1)
  initial_lambda <- NULL
  if (dist == "loglogistic"){
    initial_lambda <- log(sqrt(var(Y)*3/(pi^2)))
    selected_llk <- logistic_llk
  }else{# weibull
    initial_lambda <- log(sqrt(var(Y)*6/(pi^2)))
    selected_llk <- gumbel_llk
  }
  initial_param <- c(initial_intercept, initial_beta, initial_lambda)
  num_params <- length(initial_param)
  estim_params <- NULL
  if (!Firth || ncol(X) - length(fixed) == 1){ ## no need for Firth penalty
    estim_result <- optim(par=initial_param, fn=selected_llk,
                          Y=Y, Delta=Delta, X=X, Firth=F, fixed=fixed,
                          control=list(fnscale=-1), method="BFGS")
    estim_params <- estim_result$par
  } else{
    ## first estimate everything with Firth penalty
    estim_result1 <- optim(par=initial_param, fn=selected_llk,
                           Y=Y, Delta=Delta, X=X, Firth=T, fixed=fixed,
                           control=list(fnscale=-1), method="BFGS")
    estim_params <- estim_result1$par

    ## re-estimate the intercept and scale parameter without Firth penalty
    coefs <- estim_params[2:(length(estim_params)-1)]
    offsets <- as.vector(X[, -1, drop=F] %*% coefs)
    new_Y <- Y - offsets
    estim_result2 <- optim(par=initial_param, fn=selected_llk,
                           Y=new_Y, Delta=Delta, X=X, Firth=F, fixed=seq(2, num_params-1),
                           control=list(fnscale=-1), method="BFGS")
    estim_params[1] <- estim_result2$par[1]
    estim_params[num_params] <- estim_result2$par[num_params]
  }
  return(estim_params)
}


# LRT test
LRT_censored_regression <- function(Y, Delta, X, dist=c("loglogistic", "weibull"), Firth=T, test_param=2){

  dist <- match.arg(dist)

  full_estim_param <- NULL
  reduced_estim_param <- NULL
  teststat <- NULL
  pval <- NULL

  selected_llk <- NULL
  if (dist == "loglogistic"){
    selected_llk <- logistic_llk
  }else{# weibull
    selected_llk <- gumbel_llk
  }


  full_estim_param <- estim_censored_regression(Y=Y, Delta=Delta, X=X, Firth=Firth,
                                                 fixed=NULL, dist=dist)
  full_llk <- selected_llk(theta=full_estim_param,
                           Y=Y, Delta=Delta, X=X, Firth=F)
  reduced_estim_param <- estim_censored_regression(Y=Y, Delta=Delta, X=X, Firth=Firth,
                                                   fixed=test_param, dist=dist)
  reduced_llk <- selected_llk(theta=reduced_estim_param,
                              Y=Y, Delta=Delta, X=X, Firth=F)

  teststat <- 2*(full_llk - reduced_llk)
  pval <- 1-pchisq(teststat, df=length(test_param))


  names(full_estim_param)[1:ncol(X)] <- colnames(X)
  names(full_estim_param)[ncol(X)+1] <- "log_Scale"
  names(reduced_estim_param)[1:ncol(X)] <- colnames(X)
  names(reduced_estim_param)[ncol(X)+1] <- "log_Scale"
  result <- list(full_estim_param=full_estim_param,
                 reduced_estim_param=reduced_estim_param,
                 teststat=teststat,
                 pval=pval)

  return(result)
}


#----------- the rest are the testing code -----------------#

# library(survival)
#
# #test data
# set.seed(2023)
# time1 <- rexp(n=20, rate=0.2)
# time2 <- rexp(n=20, rate=0.1)
# alltimes <- c(time1, time2)
# event <- sample(c(0, 1), size=40, replace=TRUE, prob=c(0.7, 0.3))
# alltimes[!event] <- alltimes[!event] / 2
# neglogtime <- -log(alltimes)
# covariate <- c(rep(0, 20), rep(1, 20))
# adjust_covar <- rnorm(n=40)
# Xmat <- cbind(1, covariate, adjust_covar)
#
#
#
#
# # AFT model with log logistic regression
# loglogistic_vanilla_result <- LRT_censored_regression(Y=neglogtime, Delta=event, X=Xmat, dist="loglogistic",
#                                                    Firth=F)
#
# loglogistic_vanilla_hessian <- logistic_hessian(theta=loglogistic_vanilla_result$full_estim_param, Y=neglogtime, Delta=event, X=Xmat)
# covar_loglogistic_vanilla <- solve(-loglogistic_vanilla_hessian)
#
#
# ## run the standard survreg for comparison
# standard_full <- survreg(Surv(alltimes, event) ~ covariate + adjust_covar, dist="loglogistic")
# ## verify coefficient estimate
# stopifnot(abs(loglogistic_vanilla_result$full_estim_param[1:3] + standard_full$coefficients) < 1e-3)
# ## verify variance estimate
# stopifnot(abs(diag(covar_loglogistic_vanilla) - diag(standard_full$var)) < 1e-3)
# standard_reduced <- survreg(Surv(alltimes, event) ~ adjust_covar, dist="loglogistic")
# standard_LRTtest <- anova(standard_reduced, standard_full)
# ## verify p value
# stopifnot(abs(standard_LRTtest$`Pr(>Chi)`[2] - loglogistic_vanilla_result$pval) < 1e-3)
#
#
# loglogistic_firth_result <- LRT_censored_regression(Y=neglogtime, Delta=event, X=Xmat, dist="loglogistic",
#                                                     Firth=T)
#
#
# # AFT model with weibull distribution
#
# weibull_vanilla_result <- LRT_censored_regression(Y=neglogtime, Delta=event, X=Xmat, dist="weibull",
#                                                       Firth=F)
# weibull_vanilla_hessian <- gumbel_hessian(weibull_vanilla_result$full_estim_param, Y=neglogtime, Delta=event, X=Xmat)
# covar_weibull_vanilla <- solve(-weibull_vanilla_hessian)
#
# ## run the standard survreg for comparison
# standard_full <- survreg(Surv(alltimes, event) ~ covariate + adjust_covar, dist="weibull")
# ## verify coefficient estimate
# stopifnot(abs(weibull_vanilla_result$full_estim_param[1:3] + standard_full$coefficients) < 1e-3)
# ## verify variance estimate
# stopifnot(abs(diag(covar_weibull_vanilla) - diag(standard_full$var)) < 1e-3)
# standard_reduced <- survreg(Surv(alltimes, event) ~ adjust_covar, dist="weibull")
# standard_LRTtest <- anova(standard_reduced, standard_full)
# ## verify p value
# stopifnot(abs(standard_LRTtest$`Pr(>Chi)`[2] - weibull_vanilla_result$pval) < 1e-3)
#
#
# weibull_firth_result <- LRT_censored_regression(Y=neglogtime, Delta=event, X=Xmat, dist="weibull",
#                                                 Firth=T)


