library(lme4)
library(DiffRatio)

# case study

rm(list=ls())
# load data
folder = 'RMD/CAARS_data'
load(file.path(folder, 'CAARS_data_class_level.Rdata'))
source('RMD/ANCOM.R')

# remove extra spaces in the taxa names
taxa_names(phylodata_class) = gsub(" ", '', taxa_names(phylodata_class))
taxa_names(phylodata_class) = gsub("\\[|\\]", '', taxa_names(phylodata_class))
taxa_names(phylodata_class) = gsub("\\(|\\)", '', taxa_names(phylodata_class))
taxa_names(phylodata_class) = gsub("-", '_', taxa_names(phylodata_class))

# make sure that the taxa are the rows and the individuals are columns
count_table = otu_table(phylodata_class) %>% as.data.frame()
count_table = as.data.frame(count_table)
library_size <-colSums(count_table)

# set up sample id (repeated measurements for individuals)
sample_info = sample_data(phylodata_class) %>% as.data.frame()
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

taxa_w_struczero <- rowSums(struc_zero)
retained_rows <- which(taxa_w_struczero == 0)

filtered_count_no_struc_zero <- filtered_count[retained_rows, ]
filtered_count_no_struc_zero[is.na(filtered_count_no_struc_zero)] <- 0
filtered_count_no_struc_zero <- t(filtered_count_no_struc_zero) %>% as.data.frame()


write.csv(filtered_count_no_struc_zero,
          file.path(folder, 'filtered_count_class_level.csv'))
write.csv(filtered_metadata,
          file.path(folder, 'filtered_metadata_class_level.csv'), row.names=FALSE)



####---------------------------------------
# Comparison of different algorithms

rm(list=ls())
data_folder <- 'RMD/CAARS_data'
model_folder <- 'RMD/CAARS_Model_Summary'
counts_table <- read.csv(file.path(data_folder, 'filtered_count_class_level.csv'),
                         row.names = 1)
metadata <- read.csv(file.path(data_folder, 'filtered_metadata_class_level.csv'))

taxa_names <- colnames(counts_table)
taxa_pairs <- combn(taxa_names, 2) %>% t()

pvals_result <- as.data.frame(taxa_pairs)
colnames(pvals_result) <- c("Taxa1", "Taxa2")
pvals_result$wilx <- 0
pvals_result$ttest <- 0

# raise all values by one
shifted_counts_table <- counts_table + 1
asthma_info <- metadata$asthma
sample_id <- metadata$SAMPLE_ID

for (j in 1:nrow(pvals_result)){
  # generate a subset dataframe
  t1_name <- pvals_result$Taxa1[j]
  t2_name <- pvals_result$Taxa2[j]
  subset_data <- cbind(shifted_counts_table[, t1_name],
                       shifted_counts_table[, t2_name],
                       asthma_info) %>% as.data.frame()
  colnames(subset_data) <- c(t1_name, t2_name, "Asthma")
  subset_data$logratio <- log(subset_data[, t1_name])-log(subset_data[, t2_name])

  # wilcoxon test
  wilx_test <- wilcox.test(logratio~Asthma, data=subset_data)
  pvals_result$wilx[j] <- wilx_test$p.value

  # t test
  t_test <- t.test(logratio~Asthma, data=subset_data)
  pvals_result$ttest[j] <- t_test$p.value
}

pvals_result$GLMM_PIRLS <- 0
pvals_result$GLMM_LA_NM <- 0
pvals_result$GLMM_LA_BOBYQA <- 0
pvals_result$GLMM_GQ_NM <- 0
pvals_result$GLMM_GQ_BOBYQA <- 0

# keep track of taxa pairs that have convergence failure warning
PIRLS_failedpairs = list()
LA_NM_failedpairs = list()
LA_BOBYQA_failedpairs = list()
GQ_NM_failedpairs = list()
GQ_BOBYQA_failedpairs = list()

for (j in 1:nrow(pvals_result)){
  print(j)
  t1 <- pvals_result$Taxa1[j]
  t2 <- pvals_result$Taxa2[j]

  # penalized iteratively reweighted least squares
  tryCatch({
    glmm_PIRLS_model <- dfr(count_table=counts_table, sample_info=metadata,
                     covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = FALSE,
                     nAGQ = 0L) |> summary()
  },
  warning = function(cond){
    PIRLS_failedpairs[[length(PIRLS_failedpairs) + 1]] <<- c(t1, t2)
  },
  finally = {
    pvals_result$GLMM_PIRLS[j] <- glmm_PIRLS_model$coefficients[[8]]
  })

  # nelder mead, LA
  tryCatch({
    glmm_NM_LA <- dfr(count_table=counts_table, sample_info=metadata,
                     covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = FALSE,
                     nAGQ = 1L, optimizer = "Nelder_Mead") |> summary()
  },
  warning = function(cond){
    LA_NM_failedpairs[[length(LA_NM_failedpairs) + 1]] <<- c(t1, t2)
  },
  finally = {
    pvals_result$GLMM_LA_NM[j] <-glmm_NM_LA$coefficients[[8]]
  })

  # nelder mead, nagq=10
  tryCatch({
    glmm_NM_GQ <- dfr(count_table=counts_table, sample_info=metadata,
                      covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = FALSE,
                      nAGQ = 50L, optimizer = "Nelder_Mead") |> summary()
  },
  warning = function(cond){
    GQ_NM_failedpairs[[length(GQ_NM_failedpairs) + 1]] <<- c(t1, t2)
  },
  finally = {
    pvals_result$GLMM_GQ_NM[j] <- glmm_NM_GQ$coefficients[[8]]
  })

  # bobyqa, LA
  tryCatch({
    glmm_bobyqa_LA <- dfr(count_table=counts_table, sample_info=metadata,
                         covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = FALSE,
                         nAGQ = 1L, optimizer = "bobyqa") |> summary()
  },
  warning = function(cond){
    LA_BOBYQA_failedpairs[[length(LA_BOBYQA_failedpairs) + 1]] <<- c(t1, t2)
  },
  finally = {
    pvals_result$GLMM_LA_BOBYQA[j] <- glmm_bobyqa_LA$coefficients[[8]]
  })

  # bobyqa, nagq=10
  tryCatch({
    glmm_bobyqa_GQ <- dfr(count_table=counts_table, sample_info=metadata,
                          covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = FALSE,
                          nAGQ = 50L, optimizer = "bobyqa") |> summary()
  },
  warning = function(cond){
    GQ_BOBYQA_failedpairs[[length(GQ_BOBYQA_failedpairs) + 1]] <<- c(t1, t2)
  },
  finally = {
    pvals_result$GLMM_GQ_BOBYQA[j] <- glmm_bobyqa_GQ$coefficients[[8]]
  })
}

# identify problematic pairs of taxa through variance of p values
GLMM_results <- log(pvals_result[, seq(6, 9)])
GLMM_results$variance <- apply(GLMM_results, 1, var)
tpair_pvalvar_rank <- length(GLMM_results$variance) + 1 - rank(GLMM_results$variance)
rowids <- which(tpair_pvalvar_rank %in% seq(1, 10))
selected_pairs <- pvals_result[rowids, ]



write.csv(pvals_result, file.path('RMD/CAARS_Model_Summary', 'Pvals_class_level.csv'),
          row.names=FALSE)

failed_pairs <- list(LA_NM = LA_NM_failedpairs,
                     LA_BOBYQA = LA_BOBYQA_failedpairs,
                     GQ_NM = GQ_NM_failedpairs,
                     GQ_BOBYQA = GQ_BOBYQA_failedpairs)

saveRDS(failed_pairs,
        file=file.path('RMD/CAARS_Model_Summary', 'failed_pairs_order_level.rds'))


