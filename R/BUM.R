
# inspired by ClassComparison package

expit <- function(a) { exp(a)/(1+exp(a)) }

# log likelihood of mixture with parameters a and L
BUM_llk <- function(vec, pvals) {
  psi <- vec[1]
  phi <- vec[2]
  a <- expit(psi)
  L <- expit(phi)
  sum(log(L+(1-L)*a*pvals^(a-1)))
}

BUM_fit <- function(pvals, ...) {
  if (all(is.na(pvals))) {
    stop("all p-values were NA; nothing to compuite")
  }
  orig.pvals <- pvals
  if (any(is.na(pvals))) {
    pvals <- pvals[!is.na(pvals)]
  }
  if (any(pvals < 0) || any(pvals > 1)) {
    stop("all p-values must be between 0 and 1")
  }
  if(min(pvals)==0) {
    min.nonzero <- min(pvals[pvals>0])
    pvals[pvals==0] <- min.nonzero/2
  }
  fitted_values <- optim(c(1/2, 1/2), BUM_llk, ...,
                    method='BFGS', pvals=pvals, control=list(fnscale=-1))	# least squares fit
  psi <- fitted_values$par[1]
  phi <- fitted_values$par[2]
  ahat <- expit(psi)	# MLE estimate of a
  lhat <- expit(phi)	# MLE estimate of L
  pihat <- lhat + (1-lhat)*ahat	# upper bound on percent unchanged
  output <- list(estim_params = fitted_values$par,
                 pihat=pihat,
                 pvals=orig.pvals)
  return(output)
}


