
#' dfr
#'
#' Differential Ratio Analysis
#' @import stats
#' @param phyloseq_data phyloseq data
#' @param covar name of covariates of interest
#' @param adjust names of covariates of no interest
#' @param tpair names of taxa pairs
#' @return logistic regression result
#' @export
dfr = function(phyloseq_data, tpair, covar, adjust=c()){

  # extract the OTU counts and the sample information from the phyloseq data
  otu_counts <- phyloseq::otu_table(phyloseq_data)
  if (otu_counts@taxa_are_rows){# transpose the count data if the taxa are rows
    otu_counts <- t(otu_counts)
  }
  sample_info <- phyloseq::sample_data(phyloseq_data)
  data <- cbind(otu_counts[, tpair], sample_info[, c(covar, adjust)]) |> as.data.frame()
  data$sampleid <- seq(1, nrow(data)) # add a column representing sample id

  # set up the regression formula
  response <- sprintf("cbind(%s, %s) ~ ", tpair[1], tpair[2])
  covariates <- paste(c(covar, adjust), collapse='+')
  model <- NULL
  reffect <- '(1|sampleid)'
  regfml <- formula(sprintf("%s %s + %s", response, covariates, reffect))
  model <- lme4::glmer(regfml, data=data, family="binomial")
  invisible(model)
}
