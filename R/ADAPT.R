

#' Pooling Tobit Models for microbiome differential abundance analysis
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
#' @returns a DAresult type object contains the input and the output. Use summary and plot to explore the output
#' @export
adapt <- function(input_data, cond.var, base.cond = NULL, adj.var=NULL, censor=1,
                 prev.filter=0.05, depth.filter=1000, alpha=0.05){
  
  # preprocess input phyloseq object
  preprocessed_output <- preprocess(input_data, cond.var, base.cond, adj.var, prev.filter, depth.filter)
  
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
      distance2med <- abs(estimated_effect - median(estimated_effect, na.rm=T))
      sorted_distance <- sort(distance2med)
      ordered_taxanames <- names(sorted_distance)
      reftaxa <- ordered_taxanames[1:(length(ordered_taxanames)/2)]
    } else{
      break
    }
  }
  cat(sprintf("%d taxa selected as reference\n", length(reftaxa)))
  all_CR_results <- count_ratio(count_table=count_table, design_matrix=complete_design_matrix,
                                censor = censor, reftaxa=reftaxa, test_all=T)

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

