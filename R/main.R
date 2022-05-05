
#' dfr
#'
#' Differential Ratio Analysis
#' @import stats
#' @param count_table taxa count matrix
#' @param sample_info sample information dataframe
#' @param covar name of covariates of interest
#' @param adjust names of covariates of no interest
#' @param tpair names of taxa pairs
#' @param reff random effect level
#' @param taxa_are_rows indicator of whether the count table has taxa on row
#' @return logistic regression result
#' @export
dfr = function(count_table, sample_info, tpair, covar, adjust=c(), reff=NULL, taxa_are_rows=FALSE){

  # extract the OTU counts and the sample information from the phyloseq data
  # otu_counts <- phyloseq::otu_table(phyloseq_data)
  if (taxa_are_rows){# transpose the count data if the taxa are rows
    count_table <- t(count_table)
  }
  data <- cbind(count_table[, tpair], sample_info[, c(covar, adjust)]) |> as.data.frame()
  # set up the regression formula
  response <- sprintf("cbind(%s, %s) ~ ", tpair[1], tpair[2])
  covariates <- paste(c(covar, adjust), collapse='+')
  model <- NULL
  if (is.null(reff)){
    # use sample id by default
    data$sampleid <- seq(1, nrow(data))
    reffect <- '(1|sampleid)'
  } else{
    reffect <- sprintf('(1|%s)', reff)
  }
  regfml <- formula(sprintf("%s %s + %s", response, covariates, reffect))
  model <- lme4::glmer(regfml, data=data, family="binomial")
  invisible(model)
}
