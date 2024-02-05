
rm(list=ls())
library(ADAPT)
library(dplyr)

load("R/sysdata.rda")
count_table <- simulated_null_data$otu_table
metadata <- simulated_null_data$metadata
rownames(count_table) <- sprintf("ASV%d", seq(1, 500))
colnames(count_table) <- sprintf("Indv%d", seq(1, 129))
rownames(metadata) <- sprintf("Indv%d", seq(1, 129))
null_data <- phyloseq(otu_table(count_table, taxa_are_rows = T), sample_data(metadata))
save(null_data, file="R/sysdata.rda")

ecc16s_12 <- readRDS("experiment/phyasv_visit12.rds")
ecc_metag <- readRDS("experiment/phy_metag_plaque.rds")

metadata <- sample_data(ecc_metag)

usethis::use_data(ecc16s_12)
usethis::use_data(ecc_metag)

null_results <- adapt(input_data=null_data, cond.var="X1", adj.var=c("X2"))

realdata_results_1 <- adapt(input_data=ecc16s_12, cond.var="CaseEver", adj.var=c("geo_loc_name", "host_sex"),
                 ref.cond="Control")

realdata_results_2 <- adapt(input_data=ecc_plaque_metag, cond.var="CaseStatus",
      ref.cond="Control")

