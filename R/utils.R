# class and methods for ADAPT result


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
