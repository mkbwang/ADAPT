
rm(list=ls())
library(ADAPT)
load("R/sysdata.rda")


null_otu_table <- simulated_null_data$otu_table
null_metadata <- simulated_null_data$metadata

# use of ADAPT
begin <- proc.time()
result_boot <- adapt(otu_table = null_otu_table, metadata=null_metadata,
                      covar="X1", adjust="X2", taxa_are_rows = TRUE, boot=T)
duration <- proc.time() - begin


result_noboot <- adapt(otu_table = null_otu_table, metadata=null_metadata,
                       covar="X1", adjust="X2", taxa_are_rows = TRUE, boot=F)



# scale estimation

null_otu_table <- t(null_otu_table)
Deltamat <- 1*(null_otu_table > 0)
refcounts <- rowSums(null_otu_table)
designmat <- cbind(1, null_metadata$X1, null_metadata$X2) |> as.matrix()

scale_estimates <-  boot_estim(null_otu_table, refcounts, Deltamat, designmat,
                               boot_replicate=1000, n_boot_taxa=100)

