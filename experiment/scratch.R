
rm(list=ls())
library(ADAPT)
library(dplyr)


ecc16s_12 <- readRDS("experiment/phyasv_visit12.rds")

adapt(input_data=ecc16s_12, cond.var="CaseEver", adj.var=c("geo_loc_name", "host_sex"),
                 ref.cond="Control")

ecc_plaque_metag <- readRDS("experiment/phy_metag_plaque.rds")
result <- adapt(input_data=ecc_plaque_metag, cond.var="CaseStatus",
      ref.cond="Control")

