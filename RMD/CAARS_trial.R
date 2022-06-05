
library(DiffRatio)
library(biomformat)
library(phyloseq)
library(doParallel)
library(dplyr)
library(nlme)
library(cowplot)
library(ggplot2)
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
taxa_pairs$dfr <- 0
# count_table <- t(filtered_count)
# data <- cbind(count_table[, c(t1, t2)], filtered_metadata[, 'asthma']) %>% as.data.frame()

dfr_start <- proc.time()
for (j in 1:nrow(taxa_pairs)){
  setTxtProgressBar(pb, j)
  t1 <- taxa_pairs$Taxa1[j]
  t2 <- taxa_pairs$Taxa2[j]

  glmm_result <- dfr(count_table=filtered_count, sample_info=filtered_metadata,
                       covar=c("asthma"), tpair=c(t1, t2), reff="SAMPLE_ID", taxa_are_rows = TRUE,
                     nAGQ = 10L, optimizer = "bobyqa") |> summary()

  taxa_pairs$dfr[j] <- glmm_result$coefficients[[8]]
}
dfr_end <- proc.time()

save(res, glmm_result, file="RMD/CAARS_data/models.RData")

# taxa_pairs$ANCOM <- p.adjust(taxa_pairs$ANCOM, method='BH')
taxa_pairs$dfr_adjusted <- p.adjust(taxa_pairs$dfr, method='BH')

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


