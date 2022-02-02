library(lme4)
#' dfr
#'
#' Differential Ratio Analysis
#' @import stats
#' @param data dataframe
#' @param covar name of covariates
#' @param tpair names of taxa pairs
#' @param indveff name of random effect level, optional
#' @return logistic regression result
#' @export
dfr = function(data, covar, tpair, indveff=NULL){
  response <- sprintf("cbind(%s, %s) ~ ", tpair[1], tpair[2])
  covariates <- paste(covar, collapse='+')
  model <- NULL
  if (is.character(indveff)){ # the user offered individual level random effects
    reffect <- sprintf('(1|%s)', indveff)
    regfml <- formula(sprintf("%s %s + %s", response, covariates, reffect))
    model <- lme4::glmer(regfml, data=data, family="binomial")# glmm
  } else{
    regfml <- formula(sprintf('%s%s', response, covariates))
    model <- glm(regfml, data=data, family=binomial(link = "logit")) # glm
  }
  invisible(model)
}
