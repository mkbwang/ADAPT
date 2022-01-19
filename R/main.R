
#' dfr
#'
#' Differential Ratio Analysis
#' @import stats
#' @param data matrix or dataframe
#' @param covar name of covariates
#' @param tpair names of taxa pairs
#' @return logistic regression result
#' @export
dfr = function(data, covar, tpair){
  covar_mat = cbind(1, data[,covar])
  output = data[,tpair]
  ratios = output[,1]/(output[,1]+output[,2])
  reg_result <- glm.fit(x=covar_mat, y=ratios, family=quasibinomial(link = "logit"))
  invisible(reg_result)
}
