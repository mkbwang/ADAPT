
library(foreach)
library(doParallel)
library(parallelly)
library(dplyr)



# From dispmod package, Williams etal 1982
glm.binomial.disp <- function(object, maxit = 30, verbose = TRUE)
{
  if (class(object)[1] != "glm")
    stop("first argument must be a fitted model of class \"glm\" !")
  class <- class(object)
  if (!(family(object)$family == "binomial" & family(object)$link == "logit"))
    stop("overdispersed model fitting available only for binomial regression models with logit link function!")

  pearson.X2 <- function(x) sum(residuals(x, "pearson")^2)

  y <- object$model[,1]       # observed proportion of success & failures
  trials <- apply(y, 1, sum)  # = object$prior.weights
  X <- model.matrix(object)
  p <- length(object$coefficients)
  n <- dim(X)[[1]]
  h <- lm.influence(object)$hat
  X2 <- pearson.X2(object)
  env <- parent.frame()
  assign("object", object, envir = env)

  # initial estimate of dispersion parameter
  phi <- max((X2 - (n-p)) / sum((trials-1)*(1-h)), 0)

  if(verbose)
    cat("\nBinomial overdispersed logit model fitting...\n")

  # loop until Pearson X2 approx equal to 1
  i <- 0
  converged <- TRUE
  while( abs(X2/(n-p)-1) > object$control$epsilon )
  {
    i <- i + 1
    if(i > maxit)
    { converged <- FALSE
    break() }
    else if(verbose)
    { cat("Iter. ", i, " phi:", format(phi), "\n") }

    # computes weights
    w <- 1/(1+phi*(trials-1))
    # re-fit the model using update() evaluated in original model
    assign("disp.weights", w, envir = env)
    object <- eval(expression(update(object, weights=disp.weights)),
                   envir = env)
    #
    h <- lm.influence(object)$hat
    X2 <- pearson.X2(object)
    # current estimate of dispersion parameter
    phi <- max((X2 - sum(w*(1-h))) / sum(w*(trials-1)*(1-h)), 0)
  }

  if(verbose)
  { if(converged)
  { cat("Converged after", i, "iterations. \n")
    cat("Estimated dispersion parameter:", format(phi), "\n")
    print(summary(object))
  }
    else
      warning("algoritm not converged after ", i, " iterations!")
  }

  object <- c(object, list(dispersion = phi, disp.weights = w))
  class(object) <- class
  invisible(object)
}



reference_GLM <- function(count_data, metadata, covar, reftaxa, complement=FALSE){
  alltaxa <- rownames(count_data)
  TBDtaxa <- NULL
  if (complement){
    TBDtaxa <- setdiff(alltaxa, reftaxa)
  } else{
    TBDtaxa <- reftaxa
  }
  refcounts <- colSums(count_data[reftaxa, ]) |> as.vector()
  TBD_counts <- count_data[TBDtaxa, ]
  num_taxa <- length(TBDtaxa)
  glmdisp_result <- data.frame(Taxon = TBDtaxa,
                               effect=0,
                               SE=0,
                               teststat=0,
                               pval=0)
  
  cores=parallelly::availableCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  outcome <- foreach(j=1:nrow(glmdisp_result), .combine=rbind,
                     .packages="dplyr", .inorder=FALSE,
                     .export=c("glm.binomial.disp"),
                     .errorhandling="remove") %dopar% {
                       estimation <- list(ID=j, effect=NA, SE=NA, teststat=NA, pval=NA)
                       test_taxon <- as.character(glmdisp_result$Taxon[j])
                       selected_counts <- cbind(TBD_counts[test_taxon, ], refcounts) |> as.data.frame()
                       colnames(selected_counts) <- c('Testtaxon', 'Reftaxa')
                       selected_counts[, covar] <- metadata[, covar]
                       selected_counts$ID <- row.names(selected_counts)

                       # remove the samples with zero counts for both numerator and denominator
                       selected_counts$totcounts <- selected_counts$Reftaxa + selected_counts$Testtaxon
                       selected_counts <- selected_counts %>% filter(totcounts > 0)

                       glm_formula <- sprintf("cbind(Testtaxon, Reftaxa) ~ %s", covar) |> as.formula()
                       result <- tryCatch(
                         {
                           rawglm <- glm(glm_formula, family=binomial(link = "logit"),
                                         data=selected_counts)
                           corrected_glm <- glm.binomial.disp(rawglm, verbose = FALSE) %>%
                             summary()
                           corrected_effect <- corrected_glm$coefficients[2, 1]
                           corrected_SE <- corrected_glm$coefficients[2, 2]
                           corrected_teststat <- corrected_effect / corrected_SE
                           corrected_pval <- corrected_glm$coefficients[2, 4]
                           estimation <- list(ID=j, effect=corrected_effect,
                                              SE=corrected_SE,
                                              teststat=corrected_teststat,
                                              pval=corrected_pval)

                         },
                         error=function(cond){
                           estimation <- list(ID=j, effect=NA, SE=NA, teststat=NA, pval=NA)
                         },
                         warning=function(cond){
                           estimation <- list(ID=j, effect=NA, SE=NA, teststat=NA, pval=NA)
                         }
                       )
                       return(as.data.frame(estimation))
                     }

  stopCluster(cl)
  glmdisp_result$effect[outcome$ID] <- outcome$effect
  glmdisp_result$SE[outcome$ID] <- outcome$SE
  glmdisp_result$teststat[outcome$ID] <- outcome$teststat
  glmdisp_result$pval[outcome$ID] <- outcome$pval

  return(glmdisp_result)
}
