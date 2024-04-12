# class and methods for ADAPT result

#' Differential abundance analysis result
#' 
#' @description
#' An S4 class to represent ADAPT analysis results
#' 
#' @details
#' The analysis result object contains the analysis name, reference taxa, DA taxa,
#' detailed analysis results as a dataframe and the input phyloseq object.
#' The analysis name contains the condition variable. The reference taxa `reference`
#' must be nonempty. DA taxa `signal` may be an empty string if no taxa are differentially abundant.
#' The `details` dataframe contains the taxa names, the prevalence of taxa, the estimated
#' log10 absolute abundance fold changes, the raw hypothesis test p-values and BH-adjusted p-values.
#' 
#' @slot DAAname The name of differential abundance analysis
#' @slot reference A vector of taxa names corresponding to all the reference taxa
#' @slot signal A vector of taxa names corresponding to all the DA taxa
#' @slot details A dataframe with the analysis results for all taxa
#' @slot input Input phyloseq object
#' @importFrom methods setClass
#' @export
setClass("DAresult",
         slots=c(
           DAAname="character",
           reference="character",
           signal="character",
           details="data.frame",
           input="phyloseq"
         ))

setValidity("DAresult", function(object) {
  if (!is.character(object@DAAname)) {
    "@DAAname must be a string"
  } else if (length(object@reference) == 0) {
    "There must be at least one reference taxa"
  } else{
    TRUE
  }
})



#' @importFrom phyloseq nsamples
setMethod("show", "DAresult", function(object){
  
  cat(sprintf("Differential abundance analysis for %s\n", object@DAAname))
  cat(sprintf("%d taxa and %d samples analyzed\n", nrow(object@details), nsamples(object@input)))
  num_DAtaxa <- length(object@signal)
  if (num_DAtaxa == 1 & object@signal[1] == ""){
    num_DAtaxa <- 0
  }
  cat(sprintf("%d taxa are selected as reference and %d taxa are differentially abundant\n",
              length(object@reference), num_DAtaxa))
  
})

#' Summary of differential abundance analysis 
#'
#' @description
#' Summary function for `DAresult` type object
#' 
#' @details 
#' This customized summary function reports the dimension of input count table,
#' number of reference taxa and number of differentially abundant taxa. It also returns
#' a data frame with the detailed analysis result and taxonomy of all the taxa. The
#' user can choose to only get the detailed analysis result of DA taxa or the reference taxa
#' through the `select` parameter.
#' 
#' @param object analysis result in `DAresult` type
#' @param select Taxa whose results to be returned, can be all the taxa ("all"),
#' only the differentially abundant taxa ("da") or reference taxa ("ref").
#' @returns A dataframe with detailed analysis results
#' @export
#' @importFrom phyloseq nsamples tax_table
#' @examples
#' data(ecc_saliva)
#' saliva_results <- adapt(input_data=ecc_saliva, cond.var="CaseStatus", 
#'        base.cond="Control", adj.var="Site")
#' summary(saliva_results, select="da")
setMethod("summary", "DAresult", function(object, select=c("all", "da", "ref")){
  
  select <- match.arg(select)
  # first show overall information
  show(object)
  
  # load all the taxonomical information from the input data if any
  full_taxonomy_information = tryCatch({
    as.data.frame(tax_table(object@input))
  }, warning = function(w){
    NULL
  }, error = function(e){
    NULL
  }, finally = function(f){
    NULL
  })
  
  taxonomy_information <- NULL
  subset_result <- NULL
  combined_result <- NULL
  
  # load the full analysis result
  full_result <- object@details
  DAtaxa <- object@signal
  reftaxa <- object@reference
  if (select == "da" & length(DAtaxa) == 1 & DAtaxa[1] == ""){
    combined_result <- NULL
  } else if(select == "da"){
    subset_result <- full_result[DAtaxa, , drop=FALSE]
    if (!is.null(full_taxonomy_information)) {
      taxonomy_information <- full_taxonomy_information[DAtaxa, , drop=FALSE]
    }
  } else if(select == "ref") {
    subset_result <- full_result[reftaxa, ,drop=FALSE]
    if (!is.null(full_taxonomy_information)) {
      taxonomy_information <- full_taxonomy_information[reftaxa, , drop=FALSE]
    }
  } else{ # select all the analyzed taxa by default
    subset_result <- full_result
    if (!is.null(full_taxonomy_information)) {
      taxonomy_information <- full_taxonomy_information[rownames(subset_result), , drop=FALSE]
    }
  }
  
  # column combine the analysis result and taxonomy information
  if (!is.null(taxonomy_information) & !is.null(subset_result)){
    combined_result <- cbind(subset_result, taxonomy_information)
  } else if (!is.null(subset_result)){
    combined_result <- subset_result
  }
  
  invisible(combined_result)
  
})

#' Plotting differential abundance analysis results
#' 
#' @description
#' Volcano plot of ADAPT results
#' 
#' @details
#' The customized plot function for `DAresult` type object generates a volcano plot with
#' the differentially abundant taxa highlighted. The users can decide how many taxa with 
#' the smallest p-values are labeled on the plot.
#' 
#' @param x analysis result in `DAresult` type
#' @param n.label Number of taxa to label on the plot. Note that no taxa will be labeled if
#' no DA taxa.
#' @returns A ggplot object of the volcano plot 
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab element_text .data
#' @importFrom ggplot2 theme_bw theme scale_color_manual geom_vline
#' @importFrom ggrepel geom_label_repel
#' @export
#' @examples
#' data(ecc_saliva)
#' saliva_results <- adapt(input_data=ecc_saliva, cond.var="CaseStatus", 
#'        base.cond="Control", adj.var="Site")
#' plot(saliva_results, n.label=10)
setMethod("plot", "DAresult", function(x, n.label=5){
  
  n.label <- as.integer(n.label)
  stopifnot(n.label >= 0)
  
  model_details <- x@details
  ## just in case some taxa are too rare to be estimated
  model_details <- model_details[!is.na(model_details$teststat), ]
  ## calculate negative log10 p value
  model_details$neglog10pval <- -log10(model_details$pval)
  title_x <- sprintf("Log10 Fold Change for %s", x@DAAname)
  
  generated_plot <- NULL
  # first generate basic plot
  if(x@signal[1] != ""){ # has DA taxa
    model_details$isDA <- FALSE
    model_details[x@signal, "isDA"] <- TRUE
    
    if(n.label > 0){
      n.label <- min(n.label, length(x@signal))
      index_withlabels <- which(rank(model_details$pval) <= n.label)
      model_details$label <- ""
      model_details$label[index_withlabels] <- model_details$Taxa[index_withlabels]
    }
    
    generated_plot <- ggplot(model_details, aes(x=.data$log10foldchange, y=.data$neglog10pval))+
      geom_point(alpha=0.8, aes(color=.data$isDA)) +
      xlab(title_x) + ylab("-Log10 p-value") + theme_bw() + 
      theme(legend.position="none", axis.title=element_text(size=10), 
            axis.text=element_text(size=10)) + 
      scale_color_manual(values=c("#616161", "#ff0066")) +
      geom_vline(xintercept=0, linetype="dashed", color = "blue")
    
    if (n.label > 0){ # add labels
      generated_plot <- generated_plot + 
        geom_label_repel(aes(label = .data$label),
                         size=2,
                         max.overlaps = 20,
                         box.padding   = 0.35,
                         point.padding = 0.5,
                         segment.color = 'grey50')
    }
    
  } else{ # no DA taxa
    generated_plot <- ggplot(model_details, aes(x=.data$log10foldchange, y=.data$neglog10pval)) + 
      geom_point(alpha=0.8, color="#616161") + 
      xlab(title_x) + ylab("-Log10 p-value") + theme_bw() + 
      theme(axis.title=element_text(size=10), axis.text=element_text(size=10)) + 
      geom_vline(xintercept=0, linetype="dashed", color = "blue")
  }
  
  return(generated_plot)
  
})



