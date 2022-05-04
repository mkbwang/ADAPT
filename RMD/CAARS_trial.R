
library(DiffRatio)
library(biomformat)
library(phyloseq)
library(doParallel)
library(dplyr)
library(nlme)
library(igraph)

rm(list=ls())

# load data
folder = 'RMD/CAARS_data'
load(file.path(file.path(folder, 'CAARS_processed_GENUS.Rdata')))
source('RMD/ANCOM.R')

# remove extra spaces in the taxa names
taxa_names(CAARS.data.genus) = gsub(" ", '', taxa_names(CAARS.data.genus))
taxa_names(CAARS.data.genus) = gsub("\\[|\\]", '', taxa_names(CAARS.data.genus))
taxa_names(CAARS.data.genus) = gsub("\\(|\\)", '', taxa_names(CAARS.data.genus))
taxa_names(CAARS.data.genus) = gsub("-", '_', taxa_names(CAARS.data.genus))

# make sure that the taxa are the rows and the individuals are columns
count_table = otu_table(CAARS.data.genus) %>% as.data.frame() %>% t()
count_table = as.data.frame(count_table)
library_size <-colSums(count_table)

# set up sample id (repeated measurements for individuals)
sample_info = sample_data(CAARS.data.genus) %>% as.data.frame()
sample_info$indv <- ""
for (i in 1:nrow(sample_info)){
  split_result <- strsplit(sample_info$SAMPLE_ID[i], 'V') %>% unlist()
  sample_info$indv[i] <- split_result[1]
}
table(sample_info$asthma)
table(sample_info$indv)

colnames(count_table) <- sample_info$SAMPLE_ID

# filter taxa and samples based on zero values
feature_summary <-feature_table_pre_process(count_table, sample_info, sample_var="SAMPLE_ID",
                                            group_var='asthma', out_cut=0.05, zero_cut=0.9,
                                            lib_cut=5000, neg_lb=TRUE)


filtered_count <- feature_summary$feature_table
filtered_metadata <- feature_summary$meta_data
struc_zero <- feature_summary$structure_zeros

main_var = "asthma"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = "~ 1 | indv"
lme_control = list(maxIter = 100, msMaxIter = 100, opt = "optim")

# run ANCOM
res = ANCOM(filtered_count, filtered_metadata, struc_zero, main_var, p_adj_method,
            alpha, adj_formula, rand_formula, lme_control)

adjusted_pvals <- res$q_data

