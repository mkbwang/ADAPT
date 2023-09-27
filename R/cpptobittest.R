
#' @useDynLib PTDA
#' @import nloptr
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppParallel
#' @export
cpptobittest <- function(ngene=120, nsample=50, n_boot=500, prev=0.5, nthread=8, seed=10){

  set.seed(seed)
  Ymat <- matrix(0, nsample, ngene)
  observed <- matrix(0, nsample, ngene)
  for (j in 1:ngene){
    Ymat[, j] <- rnorm(n=nsample, mean=3, sd=2)
    observed[, j] <- rbinom(n=nsample, size=1, prob = prev)
  }
  Ymat[observed == 0] <- Ymat[observed == 0]  + 0.5
  X1 <- c(rep(0, nsample/2), rep(1, nsample/2))
  X2 <- rnorm(nsample)
  Xmat <- cbind(1, X1, X2)
  setThreadOptions(numThreads = nthread)
  ptm <- proc.time()
  inference_result <- cr_lrt(Ymat, observed, Xmat, ngene, nsample, n_boot)
  duration <- proc.time() - ptm

  return(list(Ymat=Ymat, Xmat=Xmat, Dmat=observed,
              result=inference_result, duration=duration))
}

