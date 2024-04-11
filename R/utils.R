# class and methods for ADAPT result

#' @importFrom methods setClass
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



setGeneric("name", function(x) standardGeneric("name"))
setGeneric("name<-", function(x, value) standardGeneric("name<-"))

setMethod("name", "DAresult", function(x) {
  x@DAAname
})

setMethod("name<-", "DAresult", function(x, value) {
  x@DAAname <- value
  validObject(x)
  return(x)
})

# TODO: show, summary and plot
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

#' @importFrom phyloseq nsamples tax_table
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
# 
# 
# setMethod("plot", "DAresult", function(x, n.label=5){
#   
# })
# 
# 


