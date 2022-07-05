
# load data
source('RMD/preprocess.R')
folder = 'RMD/CAARS_data'
load(file.path(file.path(folder, 'CAARS_processed_GENUS.Rdata')))
processed_data <- preprocess(CAARS.data.genus, 'SAMPLE_ID', "asthma")


filtered_count <- processed_data$feature_table
filtered_metadata <- processed_data$meta_data
struc_zero <- processed_data$structure_zeros


main_var = "asthma"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula =  NULL # "~ 1 | indv"
# assume that each sample comes from different individual for now
lme_control = list(maxIter = 100, msMaxIter = 100, opt = "optim")

# run ANCOM
res = ANCOM(filtered_count, filtered_metadata, struc_zero, main_var, p_adj_method,
            alpha, adj_formula, rand_formula) # , lme_control)

ANCOM_pvals <- res$p_data



# run differential ratio analysis

## focus on genuses with no structural zeros

included_genus <- names(which(rowSums(struc_zero) == 0))
taxa_pairs <- combn(included_genus, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("Taxa1", "Taxa2")
taxa_pairs$ANCOM <- 0

for (j in 1:nrow(taxa_pairs)){
  t1 <- taxa_pairs$Taxa1[j]
  t2 <- taxa_pairs$Taxa2[j]
  taxa_pairs$ANCOM[j] <- ANCOM_pvals[t1, t2]
}
taxa_pairs$ANCOM_adjusted <- p.adjust(taxa_pairs$ANCOM, method="BH")
pb = txtProgressBar(0, nrow(taxa_pairs), style = 3)
taxa_pairs$dfr_NM_0 <- 0
taxa_pairs$dfr_NM_1 <- 0
taxa_pairs$dfr_NM_10 <- 0
taxa_pairs$dfr_bobyqa_0 <- 0
taxa_pairs$dfr_bobyqa_1 <- 0
taxa_pairs$dfr_bobyqa_10 <- 0
# count_table <- t(filtered_count)
# data <- cbind(count_table[, c(t1, t2)], filtered_metadata[, 'asthma']) %>% as.data.frame()

NM_0_failedpairs = list()
NM_1_failedpairs = list()
NM_10_failedpairs = list()
bobyqa_0_failedpairs = list()
bobyqa_1_failedpairs = list()
bobyqa_10_failedpairs = list()

dfr_start <- proc.time()
for (j in 1:nrow(taxa_pairs)){
  setTxtProgressBar(pb, j)
  t1 <- taxa_pairs$Taxa1[j]
  t2 <- taxa_pairs$Taxa2[j]

  # nelder mead, nagq=0
  tryCatch({
    glmm_NM_0 <- dfr(count_table=filtered_count, sample_info=filtered_metadata,
                       covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = TRUE,
                       nAGQ = 0L, optimizer = "Nelder_Mead") |> summary()
  },
  warning = function(cond){
    NM_0_failedpairs[[length(NM_0_failedpairs) + 1]] <- c(t1, t2)
  },
  finally = {
    taxa_pairs$dfr_NM_0[j] <- glmm_NM_0$coefficients[[8]]
  })

  # nelder mead, nagq=1
  tryCatch({
    glmm_NM_1 <- dfr(count_table=filtered_count, sample_info=filtered_metadata,
                     covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = TRUE,
                     nAGQ = 1L, optimizer = "Nelder_Mead") |> summary()
  },
  warning = function(cond){
    NM_1_failedpairs[[length(NM_1_failedpairs) + 1]] <- c(t1, t2)
  },
  finally = {
    taxa_pairs$dfr_NM_1[j] <- glmm_NM_1$coefficients[[8]]
  })

  # nelder mead, nagq=10
  tryCatch({
    glmm_NM_10 <- dfr(count_table=filtered_count, sample_info=filtered_metadata,
                     covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = TRUE,
                     nAGQ = 10L, optimizer = "Nelder_Mead") |> summary()
  },
  warning = function(cond){
    NM_10_failedpairs[[length(NM_1_failedpairs) + 1]] <- c(t1, t2)
  },
  finally = {
    taxa_pairs$dfr_NM_10[j] <- glmm_NM_1$coefficients[[8]]
  })

  # bobyqa, nagq=0
  tryCatch({
    glmm_bobyqa_0 <- dfr(count_table=filtered_count, sample_info=filtered_metadata,
                     covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = TRUE,
                     nAGQ = 0L, optimizer = "bobyqa") |> summary()
  },
  warning = function(cond){
    bobyqa_0_failedpairs[[length(bobyqa_0_failedpairs) + 1]] <- c(t1, t2)
  },
  finally = {
    taxa_pairs$dfr_bobyqa_0[j] <- glmm_bobyqa_0$coefficients[[8]]
  })

  # bobyqa, nagq=1
  tryCatch({
    glmm_bobyqa_1 <- dfr(count_table=filtered_count, sample_info=filtered_metadata,
                         covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = TRUE,
                         nAGQ = 1L, optimizer = "bobyqa") |> summary()
  },
  warning = function(cond){
    bobyqa_1_failedpairs[[length(bobyqa_1_failedpairs) + 1]] <- c(t1, t2)
  },
  finally = {
    taxa_pairs$dfr_bobyqa_1[j] <- glmm_bobyqa_1$coefficients[[8]]
  })

  # bobyqa, nagq=10
  tryCatch({
    glmm_bobyqa_10 <- dfr(count_table=filtered_count, sample_info=filtered_metadata,
                         covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = TRUE,
                         nAGQ = 10L, optimizer = "bobyqa") |> summary()
  },
  warning = function(cond){
    bobyqa_10_failedpairs[[length(bobyqa_10_failedpairs) + 1]] <- c(t1, t2)
  },
  finally = {
    taxa_pairs$dfr_bobyqa_10[j] <- glmm_bobyqa_10$coefficients[[8]]
  })
}
dfr_end <- proc.time()

# identify three example pair of taxa whose p value is very different between different algorithms
pvals_mat <- taxa_pairs[, c("dfr_NM_0", "dfr_NM_1", "dfr_NM_10", "dfr_bobyqa_0",
                            "dfr_bobyqa_1", "dfr_bobyqa_10")]

pval_var <- apply(pvals_mat, 1, var)
tpair_pval_rank <- length(pval_var) + 1 - rank(pval_var)
rowids <- which(tpair_pval_rank %in% seq(1, 10))
selected_pairs <- taxa_pairs[rowids, ]

example1 <- selected_pairs[5, c(1,2)] %>% as.vector() %>% unname() %>% unlist()
example2 <- selected_pairs[6, c(1,2)] %>% as.vector() %>% unname() %>% unlist()
example3 <- selected_pairs[7, c(1,2)] %>% as.vector() %>% unname() %>% unlist()

df_example1 <- cbind(filtered_metadata[,c('SAMPLE_ID', 'asthma')], t(filtered_count[example1, ])) %>%
  as.data.frame()
colnames(df_example1) <- c("ID", "X", "Taxa1", "Taxa2")
df_example2 <- cbind(filtered_metadata[,c('SAMPLE_ID', 'asthma')], t(filtered_count[example2, ])) %>%
  as.data.frame()
colnames(df_example2) <- c("ID", "X", "Taxa3", "Taxa4")
df_example3 <- cbind(filtered_metadata[,c('SAMPLE_ID', 'asthma')], t(filtered_count[example3, ])) %>%
  as.data.frame()
colnames(df_example3) <- c("ID", "X", "Taxa5", "Taxa6")

write.csv(df_example1, 'RMD/CAARS_Model_Summary/example1.csv', row.names=FALSE)
write.csv(df_example2, 'RMD/CAARS_Model_Summary/example2.csv', row.names=FALSE)
write.csv(df_example3, 'RMD/CAARS_Model_Summary/example3.csv', row.names=FALSE)

# taxa_pairs$ANCOM <- p.adjust(taxa_pairs$ANCOM, method='BH')
taxa_pairs$dfr_adjusted <- p.adjust(taxa_pairs$dfr_NM_10, method='BH')

write.csv(taxa_pairs, 'RMD/CAARS_Model_Summary/Pvals.csv', row.names = FALSE)

ancom_positive_glmm_negative <- taxa_pairs %>% filter(ANCOM_adjusted < 0.05) %>% filter(dfr_adjusted > 0.05)
ancom_negative_glmm_positive <- taxa_pairs %>% filter(dfr_adjusted < 0.05) %>% filter(ANCOM_adjusted > 0.05)
choices <- c("No", "Yes")
asthma_info <- choices[filtered_metadata$asthma +1]
group.colors <- c(No="#177BB6", Yes="#B63817")



genplot <- function(tpairs){
  # retrieve the counts
  counts <- filtered_count[tpairs, ] %>% t()

  # calculate the log ratios
  logratio <- log(counts[, 1] + 1) - log(counts[, 2] + 1)

  # simplify the taxon names
  taxa1 <- strsplit(tpairs[1], split='_') %>% unlist() %>% tail(n=1)
  taxa2 <- strsplit(tpairs[2], split='_') %>% unlist() %>% tail(n=1)

  # set up dataframe
  information <- cbind(counts, logratio) %>% as.data.frame()
  colnames(information) <- c(taxa1, taxa2, "LogCountDiff")
  information$Asthma <- asthma_info # add asthma information
  information$Asthma <- as.factor(information$Asthma)

  scatter_plot <- ggplot(information, aes_string(x=taxa1, y=taxa2, color="Asthma")) +
    geom_point(position="jitter", size=0.7) + scale_colour_manual(values = group.colors) +
    xlab(taxa1) + ylab(taxa2)

  violin_plot <- ggplot(information, aes(x=Asthma, y=LogCountDiff)) + geom_violin() +
    stat_summary(fun=median, geom="point", size=2, color="red") +
    ylab("Log Count Difference")

  combined_plot <- plot_grid(scatter_plot, violin_plot, ncol=1)

  return(list(data = information,
    filename = sprintf('%s_%s.pdf', taxa1, taxa2),
              plotobj = combined_plot))
}



for (i in 1:nrow(ancom_positive_glmm_negative)){

  selected_tpairs <- c(ancom_positive_glmm_negative$Taxa1[i],
                       ancom_positive_glmm_negative$Taxa2[i])

  newplot <- genplot(selected_tpairs)
  file_fullname <- file.path('RMD', 'CAARS_Model_Summary',
                             'ancom_positive_glmm_negative',
                             newplot$filename)

  ggsave(file_fullname, newplot$plotobj, width=12, height=16,
         units='cm')

}

for (i in 1:nrow(ancom_negative_glmm_positive)){

  selected_tpairs <- c(ancom_negative_glmm_positive$Taxa1[i],
                       ancom_negative_glmm_positive$Taxa2[i])

  newplot <- genplot(selected_tpairs)
  file_fullname <- file.path('RMD', 'CAARS_Model_Summary',
                             'ancom_negative_glmm_positive',
                             newplot$filename)

  ggsave(file_fullname, newplot$plotobj, width=12, height=16,
         units='cm')

}


