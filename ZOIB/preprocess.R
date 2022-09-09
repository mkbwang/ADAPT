library(DiffRatio)
library(biomformat)
library(phyloseq)
library(doParallel)
library(dplyr)
library(nlme)
library(cowplot)
library(ggplot2)
library(igraph)


preprocess <- function(phylseqdata, ID_col, group_col,
                       outlier_cutoff=0.05, zero_prop_cutoff=0.9, min_lib=5000){
  # remove extra spaces in the taxa names
  taxa_names(phylseqdata) = gsub(" ", '', taxa_names(phylseqdata))
  taxa_names(phylseqdata) = gsub("\\[|\\]", '', taxa_names(phylseqdata))
  taxa_names(phylseqdata) = gsub("\\(|\\)", '', taxa_names(phylseqdata))
  taxa_names(phylseqdata) = gsub("-", '_', taxa_names(phylseqdata))

  # make sure that the taxa are the rows and the individuals are columns
  count_table = otu_table(phylseqdata)
  if (!taxa_are_rows(count_table)){
    count_table <- t(count_table)
  }
  count_table = as.data.frame(count_table)

  sample_info = sample_data(phylseqdata) %>% as.data.frame()
  colnames(count_table) <- sample_info[[ID_col]] %>% as.vector()

  # filter taxa and samples based on zero values
  feature_summary <-feature_table_pre_process(count_table, sample_info, sample_var=ID_col,
                                              group_var=group_col, out_cut=outlier_cutoff,
                                              zero_cut=zero_prop_cutoff,
                                              lib_cut=min_lib, neg_lb=TRUE)

  return(feature_summary)

}
