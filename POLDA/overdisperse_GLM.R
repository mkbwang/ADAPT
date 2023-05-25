
library(foreach)
library(doParallel)
library(dplyr)
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


relabd_GLM <- function(count_data, metadata, covar){
  taxa_names <- row.names(count_data)
  glmdisp_result <- data.frame(Taxon=taxa_names,
                               effect = 0,
                               SE = 0,
                               teststat = 0,
                               pval = 0)
  outcome <- foreach(j=1:length(taxa_names), .combine=rbind,
                     .packages="dplyr", .inorder=FALSE,
                     .export=c("glm.binomial.disp"),
                     .errorhandling="remove") %dopar% {
                       selected_taxon <- count_data[j, ] # count of selected taxon
                       library_size <- colSums(count_data) # summed counts of other taxa
                       counts_df <- cbind(selected_taxon, library_size) |> as.data.frame()
                       colnames(counts_df) <- c("Taxon", "Depths")
                       counts_df[, covar] <- metadata[, covar]
                       glm_formula <- sprintf("cbind(Taxon, Depths) ~ %s", covar) |> as.formula()
                       
                       rawglm <- glm(glm_formula, family=binomial(link = "logit"),
                                     data=counts_df) # raw GLM
                       corrected_glm <- glm.binomial.disp(rawglm, verbose = FALSE)
                       
                       coefficient_df <- summary(corrected_glm)
                       corrected_effect <- coefficient_df$coefficients[2, 1]
                       corrected_effect_SE <- coefficient_df$coefficients[2, 2]
                       effect_teststat <- corrected_effect / corrected_effect_SE
                       effect_pval <- coefficient_df$coefficients[2, 4]
                       estimation <- list(ID=j,
                                          effect=corrected_effect,
                                          SE=corrected_effect_SE,
                                          teststat=effect_teststat,
                                          pval=effect_pval)
                       return(as.data.frame(estimation))
                       
                     }

  glmdisp_result$effect[outcome$ID] <- outcome$effect
  glmdisp_result$SE[outcome$ID] <- outcome$SE
  glmdisp_result$teststat[outcome$ID] <- outcome$teststat
  glmdisp_result$pval[outcome$ID] <- outcome$pval
  
  return(glmdisp_result)
}


# pairwise_GLM <- function(count_data, metadata, covar){
#   # count_data is the count data matrix
#   # metadata contains the covar(variable of interest) and adjust_var(adjustment variable)
#   taxa_names <- row.names(count_data)
#   taxa_pairs <- combn(taxa_names, 2) |> t() |> as.data.frame()
#   colnames(taxa_pairs) <- c("T1", "T2")
# 
#   # dataframe for all the pairs
#   glmdisp_result <- taxa_pairs
#   glmdisp_result$effect <- 0
#   glmdisp_result$SE <- 0
#   glmdisp_result$pval <- 0
# 
#   outcome <- foreach(j=1:nrow(glmdisp_result), .combine=rbind,
#                      .packages="dplyr", .inorder=FALSE,
#                      .export=c("glm.binomial.disp"),
#                      .errorhandling="remove") %dopar% {
#                        estimation <- list(ID=j, effect=NA, SE=NA, pval=NA, Fail=TRUE)
#                        t1 <- as.character(glmdisp_result$T1[j])
#                        t2 <- as.character(glmdisp_result$T2[j])
#                        selected_counts <- count_data[c(t1, t2),] %>%
#                          t() %>% as.data.frame()
#                        colnames(selected_counts) <- c('T1', 'T2')
#                        selected_counts[, covar] <- metadata[, covar]
#                        selected_counts$ID <- row.names(selected_counts)
# 
#                        selected_counts$totcounts <- selected_counts$T1 + selected_counts$T2
#                        selected_counts <- selected_counts %>% filter(totcounts > 0)
# 
#                        glm_formula <- sprintf("cbind(T1, T2) ~ %s", covar) |> as.formula()
#                        result <- tryCatch(
#                          {
#                            rawglm <- glm(glm_formula, family=binomial(link = "logit"),
#                                          data=selected_counts)
#                            corrected_glm <- glm.binomial.disp(rawglm, verbose = FALSE) %>%
#                              summary()
#                            corrected_effect <- corrected_glm$coefficients[2, 1]
#                            corrected_SE <- corrected_glm$coefficients[2, 2]
#                            corrected_pval <- corrected_glm$coefficients[2, 4]
#                            estimation <- list(ID=j, effect=corrected_effect,
#                                               SE=corrected_SE,
#                                               pval=corrected_pval,
#                                               Fail=FALSE)
# 
#                          },
#                          error=function(cond){
#                            estimation <- list(ID=j, effect=NA, SE=NA, pval=NA, Fail=TRUE)
#                          },
#                          warning=function(cond){
#                            estimation <- list(ID=j, effect=NA, SE=NA, pval=NA, Fail=TRUE)
#                          }
#                        )
#                        return(as.data.frame(estimation))
#                      }
#   glmdisp_result$effect[outcome$ID] <- outcome$effect
#   glmdisp_result$SE[outcome$ID] <- outcome$SE
#   glmdisp_result$pval[outcome$ID] <- outcome$pval
# 
#   return(glmdisp_result)
# }



reference_GLM <- function(count_data, metadata, covar, reftaxa){
  alltaxa <- rownames(count_data)
  TBDtaxa <- setdiff(alltaxa, reftaxa)
  refcounts <- colSums(count_data[reftaxa, ]) |> as.vector()
  TBD_counts <- count_data[TBDtaxa, ]
  num_taxa <- length(TBDtaxa)
  glmdisp_result <- data.frame(Taxon = TBDtaxa,
                               effect=0,
                               SE=0,
                               teststat=0,
                               pval=0)
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

  glmdisp_result$effect[outcome$ID] <- outcome$effect
  glmdisp_result$SE[outcome$ID] <- outcome$SE
  glmdisp_result$teststat[outcome$ID] <- outcome$teststat
  glmdisp_result$pval[outcome$ID] <- outcome$pval

  return(glmdisp_result)
}
