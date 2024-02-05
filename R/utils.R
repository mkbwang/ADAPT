# class and methods for ADAPT result

#' @importFrom methods setClass
setClass("DAresult",
         slots=c(
           reference="character",
           signal="character",
           details="data.frame",
           input="phyloseq"
         ))


# TODO: print
# TODO: summary
# TODO: plot
