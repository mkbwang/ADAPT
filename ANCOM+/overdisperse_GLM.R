
library(foreach)
library(doParallel)
cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


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


pairwise_GLM <- function(count_data, metadata, covar){
  # count_data is the count data matrix
  # metadata contains the covar(variable of interest) and adjust_var(adjustment variable)
  taxa_names <- row.names(count_data)
  taxa_pairs <- combn(taxa_names, 2) |> t() |> as.data.frame()
  colnames(taxa_pairs) <- c("T1", "T2")

  # dataframe for all the pairs
  glmdisp_result <- taxa_pairs
  glmdisp_result$effect <- 0
  glmdisp_result$SE <- 0
  glmdisp_result$pval <- 0

  outcome <- foreach(j=1:nrow(glmdisp_result), .combine=rbind,
                     .packages="dplyr", .inorder=FALSE,
                     .export=c("glm.binomial.disp"),
                     .errorhandling="remove") %dopar% {
                       result <- list(ID=j, effect=NA, SE=NA, pval=NA)
                       t1 <- as.character(glmdisp_result$T1[j])
                       t2 <- as.character(glmdisp_result$T2[j])
                       selected_counts <- count_data[c(t1, t2),] %>%
                         t() %>% as.data.frame()
                       colnames(selected_counts) <- c('T1', 'T2')
                       selected_counts[, covar] <- metadata[, covar]
                       selected_counts$ID <- row.names(selected_counts)

                       selected_counts$totcounts <- selected_counts$T1 + selected_counts$T2
                       selected_counts <- selected_counts %>% filter(totcounts > 0)

                       glm_formula <- sprintf("cbind(T1, T2) ~ %s", covar) |> as.formula()
                       rawglm <- glm(glm_formula, family=binomial(link = "logit"),
                                     data=selected_counts)

                       corrected_glm <- glm.binomial.disp(rawglm, verbose = FALSE) %>%
                         summary()

                       result$effect <- corrected_glm$coefficients[2, 1]
                       result$SE <- corrected_glm$coefficients[2, 2]
                       result$pval <- corrected_glm$coefficients[2, 4]

                       return(as.data.frame(result))
                     }
  glmdisp_result$effect[outcome$ID] <- outcome$effect
  glmdisp_result$SE[outcome$ID] <- outcome$SE
  glmdisp_result$pval[outcome$ID] <- outcome$pval

  return(glmdisp_result)
}

