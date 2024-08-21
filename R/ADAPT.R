

#' ADAPT
#' 
#' @description 
#' Analysis of microbiome differential abundance by pooling tobit models
#' 
#' @details
#' ADAPT takes in a metagenomics count table as a phyloseq object. 
#' The phyloseq object needs to have metadata containing 
#' at least one variable `cond.var` representing the conditions that the user is testing on. 
#' The condition variable `cond.var` can be numeric (as a continuous variable) or character (representing categorical variable).
#' ADAPT does not support multigroup comparison yet. If there are multiple conditions,
#' the user can specify the condition to single out through `base.cond`. ADAPT then carry out DAA between
#' the selected `base.cond` and all the others. ADAPT allows adjusting for other covariates. The 
#' user can specify all the covariates to adjust for by specifying `adj.var` with a vector of variable names.
#'
#' Differential abundance analysis may be too challenging for rare taxa and samples with too low sequencing depth.
#' The users can filter out taxa whose prevalences are lower than `prev.filter` (default 0.05). The users 
#' can also filter out samples whose sequencing depths (library sizes) are smaller than `depth.filter` (default 1000).
#' 
#' One major feature of ADAPT is treating zero counts as left censored observations and use Tobit models for log count ratios.
#' The zero counts by default are left censored at one. The users can change the value to censor at through `censor`.
#' Change the cutoff of BH-adjusted p-values with `alpha` (default 0.05) for calling DA taxa.
#' 
#' The returned value of `adapt` is a customized S4 type called `DAresult`. We have developed two helper functions `summary`
#' and `plot` for this special data type.
#' 
#' 
#' @useDynLib ADAPT
#' @param input_data a phyloseq object
#' @param cond.var the variable representing the conditions to compare, a character string
#' @param base.cond the condition chosen as baseline. This is only used when the condition is categorical.
#' @param adj.var the names of the variables to be adjusted, a vector of character strings
#' @param censor the value to censor at for zero counts, default 1
#' @param prev.filter taxa whose prevalences are smaller than the cutoff will be excluded from analysis, default 0.05
#' @param depth.filter a sample would be discarded if its library size is smaller than the threshold
#' @param alpha the cutoff of the adjusted p values
#' @importFrom stats optim median p.adjust qchisq
#' @returns a `DAresult` type object contains the input and the output. Use summary and plot to explore the output
#' @export
#' 
#' @examples
#' data(ecc_plaque)
#' plaque_results <- adapt(input_data=ecc_plaque, cond.var="CaseStatus", 
#'        base.cond="case")
#' data(ecc_saliva)
#' saliva_results <- adapt(input_data=ecc_saliva, cond.var="CaseStatus", 
#'        base.cond="Control", adj.var="Site")
adapt <- function(input_data, cond.var, base.cond = NULL, adj.var=NULL, censor=1,
                 prev.filter=0.05, depth.filter=1000, alpha=0.05){
  
  # preprocess input phyloseq object
  preprocessed_output <- preprocess(input_data, cond.var, base.cond, adj.var, prev.filter, depth.filter)
  
  # check if the arguments are valid
  stopifnot("The zero counts need to be censored at a positive number!" = censor > 0)
  stopifnot("The cutoff for adjusted p-values needs to be between 0 and 0.5!" = alpha > 0 & alpha < 0.5)
  
  count_table <- preprocessed_output$count_table
  complete_design_matrix <- preprocessed_output$design_matrix
  DAAname <- preprocessed_output$DAAname
  
  taxa_names <- colnames(count_table)

  
  reftaxa <- taxa_names # initially all the taxa are reference taxa(relative abundance)
  cat("Selecting Reference Set... ")
  while(1){
    relabd_result <- count_ratio(count_table = count_table, design_matrix = complete_design_matrix,
                                  censor = censor, reftaxa = reftaxa, test_all=FALSE)
    estimated_effect <- relabd_result$log10foldchange
    pvals <- relabd_result$pval
    names(pvals) <- relabd_result$Taxa
    names(estimated_effect) <- relabd_result$Taxa
    # check distribution of p values
    bumfit <- BUM_fit(pvals)
    loglik <- BUM_llk(bumfit$estim_params, pvals[!is.na(pvals)])
    if (2*loglik > qchisq(0.95, 1)){ # need to continue shrinking reference taxa set
      distance2med <- abs(estimated_effect - median(estimated_effect, na.rm=TRUE))
      sorted_distance <- sort(distance2med)
      ordered_taxanames <- names(sorted_distance)
      reftaxa <- ordered_taxanames[seq_len(length(ordered_taxanames)/2)]
    } else{
      break
    }
  }
  cat(sprintf("%d taxa selected as reference\n", length(reftaxa)))
  all_CR_results <- count_ratio(count_table=count_table, design_matrix=complete_design_matrix,
                                censor = censor, reftaxa=reftaxa, test_all=TRUE)

  all_pvals <- all_CR_results$pval
  names(all_pvals) <- all_CR_results$Taxa
  all_adjusted_pvals <- p.adjust(all_pvals, method="BH")
  all_CR_results$adjusted_pval <- all_adjusted_pvals

  significant_pvals <- all_adjusted_pvals[all_adjusted_pvals < alpha & !is.na(all_adjusted_pvals)]
  if (length(significant_pvals) > 0){
    DiffTaxa <- names(significant_pvals)
  } else{
    DiffTaxa <- c()
  }
  cat(sprintf("%d differentially abundant taxa detected\n", length(DiffTaxa)))
  if(is.null(DiffTaxa)) DiffTaxa <- ""
  output <- new("DAresult",
                DAAname=DAAname,
                reference=reftaxa, 
                signal=DiffTaxa,
                details=all_CR_results,
                input=input_data)

  invisible(output)
}

