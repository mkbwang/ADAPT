library(survival)
rm(list=ls())
#test data
set.seed(2023)
time1 <- rexp(n=20, rate=0.2)
time2 <- rexp(n=20, rate=0.1)
alltimes <- c(time1, time2)
event <- sample(c(0, 1), size=40, replace=TRUE, prob=c(0.5, 0.5))
alltimes[!event] <- alltimes[!event] / 2

covariate <- c(rep(0, 20), rep(1, 20))
Xmat <- cbind(1, covariate)



# AFT model with log logistic regression
neglogtime <- -log(alltimes)
initial_beta0 <- mean(neglogtime)
initial_beta1 <- 0
initial_lambda <- log(sqrt(var(neglogtime)*3/(pi^2)))
initial_param <- c(initial_beta0, initial_beta1, initial_lambda)

logistic_llk <- function(theta, Y, Delta, X){

  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  lambda <- theta[num_params] # scale parameter
  eta <- as.vector(X %*% beta)
  z <- (eta - Y)/exp(lambda)
  exp_z_plusone <- 1 + exp(z)
  llk <- sum(Delta * z)- sum(Delta * lambda) -
    sum((1+Delta) * log(exp_z_plusone))

  return(llk)
}

logistic_estim_result <- optim(par=initial_param, fn=logistic_llk,
                      Y=neglogtime, Delta=event, X=Xmat,
                      control=list(fnscale=-1), method="BFGS",
                      hessian=T)
logistic_estimated_parameter <- logistic_estim_result$par
logistic_estimated_covar <- solve(-logistic_estim_result$hessian)


real_result_logistic <- survreg(Surv(alltimes, event) ~ covariate, dist="loglogistic")


logistic_hessian <- function(theta, Y, Delta, X){
  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  lambda <- theta[num_params] # scale parameter
  eta <- as.vector(X %*% beta)
  z <- (eta - Y)/exp(lambda)

  output_hessian <- matrix(0, nrow=num_params, ncol=num_params)

  ## first calculate the hessian for the beta only
  weights <- -(1+Delta)/exp(2*lambda)*exp(z)/(1+exp(z))^2
  output_hessian[1:(num_params-1), 1:(num_params-1)] <- t(X) %*% diag(weights) %*% X

  ## then calculate the hessian for lambda
  output_hessian[num_params, num_params] <- -sum((exp(z) - Delta)*z/(1+exp(z))) -
    sum((1+Delta)*exp(z)*z^2/(1+exp(z))^2)

  ## finally calculate derivative over both beta and lambda
  weights <- 1/exp(lambda)*(1-(1+Delta)/(1+exp(z)) + (1+Delta)*z*exp(z)/(1+exp(z))^2)
  output_hessian[1:(num_params-1), num_params] <- output_hessian[num_params, 1:(num_params-1)] <-
    as.vector(weights %*% X)

  return(output_hessian)
}


logistic_derived_hmat <- logistic_hessian(theta=logistic_estimated_parameter,
                                 Y=neglogtime, Delta=event, X=Xmat)

stopifnot(max(abs(logistic_derived_hmat - logistic_estim_result$hessian)) < 1e-2)

# AFT result with Weibull distribution

initial_beta0 <- mean(neglogtime)
initial_beta1 <- 0
initial_lambda <- log(sqrt(var(neglogtime)*6/(pi^2)))
initial_param <- c(initial_beta0, initial_beta1, initial_lambda)


gumbel_llk <- function(theta, Y, Delta, X){

  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  lambda <- theta[num_params] # scale parameter
  eta <- as.vector(X %*% beta)
  z <- (eta - Y)/exp(lambda)

  llk <- -sum(Delta * lambda) + sum(Delta * z) - sum(exp(z))
  return(llk)

}

gumbel_estim_result <- optim(par=initial_param, fn=gumbel_llk,
                               Y=neglogtime, Delta=event, X=Xmat,
                               control=list(fnscale=-1), method="BFGS",
                               hessian=T)
gumbel_estimated_parameter <- gumbel_estim_result$par
gumbel_estimated_covar <- solve(-gumbel_estim_result$hessian)

real_result_weibull <- survreg(Surv(alltimes, event) ~ covariate, dist="weibull")


gumbel_hessian <- function(theta, Y, Delta, X){

  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
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
  weights <- exp(-lambda)*(exp(z)*z+exp(z)-Delta)
  output_hessian[num_params, 1:(num_params-1)] <- output_hessian[1:(num_params-1), num_params] <-
    as.vector(weights %*% X)

  return(output_hessian)

}

gumbel_derived_hmat <- gumbel_hessian(theta=gumbel_estimated_parameter,
                               Y=neglogtime, Delta=event,X=Xmat)

stopifnot(max(abs(gumbel_derived_hmat - gumbel_estim_result$hessian)) < 1e-2)

#TODO: Add the penalty based on Firth correction

