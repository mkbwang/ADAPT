

#' Preprocess input phyloseq data to get count table and metadata design matrix
#' @param input_data a phyloseq object
#' @param cond.var the variable representing the conditions to compare, a character string
#' @param base.cond the condition chosen as baseline. This is only used when the condition is categorical.
#' @param adj.var the names of the variables to be adjusted, a vector of character strings
#' @param prev.filter taxa whose prevalences are smaller than the cutoff will be excluded from analysis
#' @param depth.filter a sample would be discarded if its library size is smaller than the threshold
#' @importFrom phyloseq filter_taxa prune_samples otu_table sample_data taxa_are_rows sample_sums
#' @importFrom stats model.matrix
#' @returns a list with count table, metadata design matrix and DAA name
preprocess <- function(input_data, cond.var, base.cond, adj.var, prev.filter, depth.filter){
  
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
    if (is.null(base.cond)){
      base.cond <- unique(main_variable)[1]
    } else if (!base.cond %in% main_variable) {
      base.cond <- unique(main_variable)[1]
    }
    cat(sprintf("Choose '%s' as the baseline condition\n", base.cond))
    if (length(unique(main_variable)) == 2){
      others <- setdiff(main_variable, base.cond)
    } else{
      others <- "others"
    }
    main_variable <- as.integer(main_variable != base.cond)
    cond.var <- sprintf("%s (%s VS %s)", cond.var, others, base.cond)
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
  
  output <- list(count_table=count_table, design_matrix=complete_design_matrix,
                 DAAname = cond.var)
  
  return(output)
  
}
