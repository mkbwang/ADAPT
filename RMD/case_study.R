library(lme4)

# case study

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

# select a specific pair

Lentimicrobium_Filifactor_count <- filtered_count[c('g__Lentimicrobium',
                                                    'g__Filifactor'), ] %>% t() %>%
  as.data.frame()
Lentimicrobium_Filifactor = cbind(Lentimicrobium_Filifactor_count, filtered_metadata)
write.csv(Lentimicrobium_Filifactor, 'RMD/CAARS_data/Lentimicrobium_Filifactor.csv',
          row.names=FALSE)

data <- read.csv('RMD/CAARS_data/Lentimicrobium_Filifactor.csv')
data_replicate <- data


## for ancom, we need to impute zeros first
lenti <- data$g__Lentimicrobium
lenti[lenti==0] <- 1
filifactor <- data$g__Filifactor
filifactor[filifactor==0] <- 1
data_replicate$g__Lentimicrobium <- lenti
data_replicate$g__Filifactor <- filifactor

data_replicate$logratio <- log(data_replicate$g__Lentimicrobium) - log(data_replicate$g__Filifactor)

## Wilcoxon rank-sum test
wilx_test <- wilcox.test(logratio~asthma, data=data_replicate)
wilx_pval <- wilx_test$p.value

## T test
t_test <- t.test(logratio ~ asthma, data=data_replicate)
t_pval <- t_test$p.value

## Logistic Regression
glm_result <- glm(cbind(g__Lentimicrobium, g__Filifactor) ~ asthma, data=data,
                  family = binomial(link = "logit"))
### hand calculation
asthma_patients <- data %>% filter(asthma == 1)
healthy_patients <- data %>% filter(asthma == 0)

healthy_lenti_count <- sum(healthy_patients$g__Lentimicrobium)
healthy_filifactor_count <- sum(healthy_patients$g__Filifactor)
asthma_lenti_count <- sum(asthma_patients$g__Lentimicrobium)
asthma_filifactor_count <- sum(asthma_patients$g__Filifactor)

glm_intercept <- log(healthy_lenti_count/
                       healthy_filifactor_count)
sd_glm_intercept <- sqrt(1/healthy_lenti_count + 1/healthy_filifactor_count)

glm_asthma_effect <- log(asthma_lenti_count/
                           asthma_filifactor_count) -
  glm_intercept
sd_glm_asthma <- sqrt(1/healthy_lenti_count + 1/healthy_filifactor_count+
                        1/asthma_lenti_count + 1/asthma_filifactor_count)


## GLMM
glmm_result <- glmer(cbind(g__Lentimicrobium, g__Filifactor) ~ asthma + (1|SAMPLE_ID), data=data,
                     family = binomial(link = "logit"), nAGQ=5L,
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e5))) # nAGQ is a big issue

## other methods: nloptwrap
