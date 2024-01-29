

#' Pooling Tobit Models for microbiome differential abundance analysis
#' @useDynLib ADAPT
#' @param input_data a phyloseq object
#' @param cond.var the variable representing the conditions to compare, a character string
#' @param ref.cond the condition chosen as baseline. This is only used when the condition is categorical.
#' @param adj.var the names of the variables to be adjusted, a vector of character strings
#' @param prev.filter taxa whose prevalences are smaller than the cutoff will be excluded from analysis, default 0.05
#' @param depth.filter a sample would be discarded if its library size is smaller than the threshold
#' @param alpha the cutoff of the adjusted p values
#' @importFrom stats optim
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats qchisq
#' @returns a DAresult type object contains the input and the output. Use summary and plot to explore the output
#' @export
adapt <- function(input_data, metadata, cond.var, ref.cond = NULL, adj.var=NULL,
                 prev.filter=0.05, depth.filter=1000, alpha=0.05){
  
  # check if input data type is phyloseq
  if (class(input_data)[1] != "phyloseq"){
    stop("Input data isn't a phyloseq object!")
  }
  # filter phyloseq object based on taxa prevalence and sequencing depth
  subset_data <- filter_taxa(input_data, function(x) mean(x>0) > prev.filter, TRUE)
  subset_data <- prune_samples(sample_sums(subset_data) > depth.filter, subset_data)
  
  # check if the variables exist in the metadata
  metadata <- data.frame(sample_data(subset_data))
  allcols <- colnames(metadata)
  if (class(cond.var) != "character"){
    stop("The main variable name for conditions need to be a string.")
  }
  if (class(adj.var) != "character" & class(adj.var) != "NULL"){
    stop("The variables for adjustments should be either NULL or a vector of character strings.")
  }
  selected_cols <- c(cond.var, adj.var)
  if (!all(selected_cols %in% allcols)){
    unavailable_cols <- selected_cols[!selected_cols %in% allcols]
    stop(sprintf("Some columns are not available in the metadata! (%s)", 
                 paste(unavailable_cols, collapse=",")))
  }
  # dichotomize categorical variables if selected, set up design matrix
  if (any(is.na(metadata))){
    stop("No missing data allowed in the metadata!")
  }
  main_variable <- metadata[, cond.var]
  if (length(unique(main_variable)) == 1){
    stop("All samples share the same condition!")
  }
  if (!is.numeric(main_variable)){
    if (is.null(ref.cond) | !ref.cond %in% main_variable){
      ref.cond <- unique(main_variable)[1]
    }
    cat(sprintf("Choose '%s' as the reference condition\n", ref.cond))
    main_variable <- main_variable != ref.cond
    if (length(unique(main_variable)) == 2){
      others <- unique(main_variable)[2]
    } else{
      others <- "others"
    }
    cond.var <- sprintf("%s_%sVS%s", cond.var, others, ref.cond)
  }
  adjustments <- NULL
  if (!is.null(adj.var)){
    adjustments <- metadata[, adj.var, drop=F]
    adjustments<- model.matrix(~., data=adjustments)
    adjustments<- adjustments[, -1] # remove intercept
  }
  complete_design_matrix <- cbind(1, main_variable, adjustments)
  
  # parse the count matrix and edit the count_ratio function
  count_table <- otu_table(subset_data)
  if (taxa_are_rows(subset_data)){
    count_table <- t(count_table)
  }

  if (any(is.na(count_table))){
    stop("No missing data allowed in the count table!")
  }
  
  cat(sprintf("%d taxa and %d samples being analyzed...\n", 
              ncol(count_table), nrow(count_table)))
  
  taxonomies <- tax_table(subset_data)
  
  taxa_names <- rownames(taxonomies)
  reftaxa <- taxa_names # initially all the taxa are reference taxa(relative abundance)
  cat("Selecting Reference Set...")
  while(1){
    relabd_result <- count_ratio(count_table = count_table, design_matrix = complete_design_matrix,
                                  reftaxa = reftaxa, test_all=FALSE)
    estimated_effect <- relabd_result$effect
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
                                reftaxa=reftaxa, test_all=T)

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
  cat(sprintf("%d DA taxa detected\n", length(DiffTaxa)))
  output <- new("DAresult", 
                reference=reftaxa, 
                signal=DiffTaxa,
                details=all_CR_results,
                input=input_data)

  invisible(output)
}

