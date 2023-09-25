
#' @useDynLib PTDA
#' @import nloptr
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppParallel
#' @export
cpptobittest <- function(){
  ngene <- 200
  nsample <- 50

  Ymat <- matrix(0, nsample, ngene)
  observed <- matrix(0, nsample, ngene)
  for (j in 1:ngene){
    Ymat[, j] <- rnorm(n=nsample, mean=3, sd=1)
    observed[, j] <- rbinom(n=nsample, size=1, prob=0.6)
  }
  Ymat[observed == 0] <- Ymat[observed == 0]  - 0.5
  X1 <- c(rep(0, nsample/2), rep(1, nsample/2))
  X2 <- rnorm(nsample)
  Xmat <- cbind(1, X1, X2)

  inference_result <- cr_lrt(Ymat, observed, Xmat, ngene, nsample)

  return(inference_result)
}

