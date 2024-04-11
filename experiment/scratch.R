
rm(list=ls())
# library(ADAPT)
library(dplyr)

load("R/sysdata.rda")

# null test case
null_output <- adapt(input_data=phyobj_null, 
                      cond.var="main", adj.var=c("confounder"))
null_results <- summary(null_output)
FPR <- mean(null_results$pval < 0.05)



# balanced test case
phyobj_balanced <- balanced_case$phyobj
balanced_truth <- balanced_case$truth
balanced_output <- adapt(input_data=phyobj_balanced, 
                      cond.var="main", adj.var=c("confounder"))
balanced_results <- summary(balanced_output, 'da')
signals <- balanced_results$Taxa
true_signals <- balanced_truth$taxaname[balanced_truth$isDA]
FDR_balanced <- 1-mean(signals %in% true_signals)



# unbalanced test case
phyobj_unbalanced <- unbalanced_case$phyobj
unbalanced_truth <- unbalanced_case$truth
unbalanced_output <- adapt(input_data=phyobj_unbalanced, 
                          cond.var="main", adj.var="confounder")
unbalanced_result <- summary(unbalanced_output, 'da')
signals <- unbalanced_result$Taxa
true_signals <- unbalanced_truth$taxaname[unbalanced_truth$isDA]
FDR_unbalanced <- 1-mean(signals %in% true_signals)



# saliva DAA
load("data/ecc_saliva.rda")
saliva_output <- adapt(input_data=ecc_saliva, cond.var="CaseStatus", censor=1)
saliva_result <- summary(saliva_output)


# plaque DAA
load("data/ecc_plaque.rda")
plaque_output <- adapt(input_data=ecc_plaque, cond.var="CaseStatus", censor=1)
plaque_result <- summary(plaque_output)


