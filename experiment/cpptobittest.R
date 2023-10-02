

#' @export
testcasesim <- function(ngene=120, nsample=50, prev=0.5, seed=10){
  set.seed(seed)
  Ymat <- matrix(0, nsample, ngene)
  observed <- matrix(0, nsample, ngene)
  for (j in 1:ngene){
    Ymat[, j] <- rnorm(n=nsample, mean=3, sd=2)
    observed[, j] <- c(rbinom(n=nsample/2, size=1, prob = prev), rbinom(n=nsample/2, size=1, prob = prev))
  }
  Ymat[observed == 0] <- Ymat[observed == 0]  + 1
  X1 <- c(rep(0, nsample/2), rep(1, nsample/2))
  X2 <- rnorm(nsample)
  X3 <- rnorm(nsample)
  Xmat <- cbind(1, X1, X2, X3)
  return(list(Ymat=Ymat, Delta=observed, Xmat=Xmat))
}


#' @useDynLib PTDA
#' @import nloptr
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppParallel
#' @export
cpptobittest <- function(simdata, n_boot=500, nthread=8){

  Ymat <- simdata$Ymat

  setThreadOptions(numThreads = nthread)

  ptm <- proc.time()
  inference_result <- cr_estim(simdata$Ymat, simdata$Delta, simdata$Xmat)
  inference_duration <- proc.time() - ptm

  ptm <- proc.time()
  boot_result <- boot_estim(Ymat, observed, Xmat, boot_replicate=n_boot)
  boot_duration <- proc.time() - ptm

  return(list(inference_result=inference_result, inference_duration=inference_duration,
              boot_result=boot_result, boot_duration=boot_duration))
}

