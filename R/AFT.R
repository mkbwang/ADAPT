library(survival)

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
initial_s <- sqrt(var(neglogtime)*3/(pi^2))
initial_param <- c(initial_beta0, initial_beta1, log(initial_s))

logistic_llk <- function(theta, Y, Delta, X){

  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  s <- exp(theta[num_params]) # scale parameter
  eta <- as.vector(X %*% beta)
  neg_resid_y <- (eta - Y)/s
  exp_res_plusone <- 1+exp(neg_resid_y)
  llk <- sum(Delta * neg_resid_y)- sum(Delta * log(s)) -
    sum((1+Delta) * log(exp_res_plusone))

  return(llk)
}

logistic_estim_result <- optim(par=initial_param, fn=logistic_llk,
                      Y=neglogtime, Delta=event, X=Xmat,
                      control=list(fnscale=-1), method="BFGS",
                      hessian=T)
logistic_estimated_parameter <- logistic_estim_result$par
logistic_estimated_covar <- solve(-logistic_estim_result$hessian)


real_result_logistic <- survreg(Surv(alltimes, event) ~ covariate, dist="loglogistic")


# AFT result with Weibull distribution

initial_beta0 <- mean(neglogtime)
initial_beta1 <- 0
initial_s <- sqrt(var(neglogtime)*6/(pi^2))
initial_param <- c(initial_beta0, initial_beta1, log(initial_s))


gumbel_llk <- function(theta, Y, Delta, X){
  stopifnot(ncol(X)==length(theta)-1) # check dimension
  num_params <- length(theta)
  beta <- theta[1:(num_params-1)] # coefficients for location
  s <- exp(theta[num_params]) # scale parameter
  eta <- as.vector(X %*% beta)
  neg_resid_y <- (eta - Y)/s
  exp_res_plusone <- 1+exp(neg_resid_y)
  llk <- -sum(Delta * log(s)) + sum(Delta * neg_resid_y) - sum(exp(neg_resid_y))
  return(llk)
}

gumbel_estim_result <- optim(par=initial_param, fn=gumbel_llk,
                               Y=neglogtime, Delta=event, X=Xmat,
                               control=list(fnscale=-1), method="BFGS",
                               hessian=T)
gumbel_estimated_parameter <- gumbel_estim_result$par
gumbel_estimated_covar <- solve(-gumbel_estim_result$hessian)

real_result_weibull <- survreg(Surv(alltimes, event) ~ covariate, dist="weibull")


#TODO: Add the penalty based on Firth correction

