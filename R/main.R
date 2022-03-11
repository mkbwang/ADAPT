library(lme4)
#' dfr
#'
#' Differential Ratio Analysis
#' @import stats
#' @param data dataframe
#' @param covar name of covariates
#' @param tpair names of taxa pairs
#' @return logistic regression result
#' @export
dfr = function(data, covar, tpair){
  data$sampleid <- seq(1, nrow(data)) # add a column representing sample id
  response <- sprintf("cbind(%s, %s) ~ ", tpair[1], tpair[2])
  covariates <- paste(covar, collapse='+')
  model <- NULL
  reffect <- '(1|sampleid)'
  regfml <- formula(sprintf("%s %s + %s", response, covariates, reffect))
  model <- lme4::glmer(regfml, data=data, family="binomial")
  invisible(model)
}
