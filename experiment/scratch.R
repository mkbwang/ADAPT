
rm(list=ls())
# library(ADAPT)
library(dplyr)

load("R/sysdata.rda")

# null test case
null_results <- adapt(input_data=phyobj_null, 
                      cond.var="main", adj.var=c("confounder"))
result1 <- null_results@details
result2 <- null_results@details
FPR <- mean(null_results@details$pval < 0.05)

# balanced test case
phyobj_balanced <- balanced_case$phyobj
balanced_truth <- balanced_case$truth
balanced_results <- adapt(input_data=phyobj_balanced, 
                      cond.var="main", adj.var=c("confounder"))
signals <- balanced_results@signal
true_signals <- balanced_truth$taxaname[balanced_truth$isDA]
FDR_balanced <- 1-mean(signals %in% true_signals)

# unbalanced test case
phyobj_unbalanced <- unbalanced_case$phyobj
unbalanced_truth <- unbalanced_case$truth
unbalanced_results <- adapt(input_data=phyobj_unbalanced, 
                          cond.var="main", adj.var="confounder")
signals <- unbalanced_results@signal
true_signals <- unbalanced_truth$taxaname[unbalanced_truth$isDA]
FDR_unbalanced <- 1-mean(signals %in% true_signals)


load("data/ecc_saliva.rda")
# saliva DAA
realdata_results_1 <- adapt(input_data=ecc_saliva, cond.var="CaseStatus", ref.cond="Control")

# plaque DAA
load("data/ecc_plaque.rda")
realdata_results_2 <- adapt(input_data=ecc_plaque, cond.var="CaseStatus", ref.cond="Case")

