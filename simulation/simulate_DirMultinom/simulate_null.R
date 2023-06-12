
remove(list=ls())

folder <- '/home/wangmk/MDAWG/POLDA/simulation'
source(file.path(folder, "simulate_DirMultinom/utils.R"))


AGP_template <- readRDS(file.path(folder, 'data', 'AGP_template.rds'))
AGP_count_mat <- otu_table(AGP_template)@.Data

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

nSample <- 100
nTaxa <- nrow(AGP_count_mat)
# permute counts to get a null dataset
AGP_permuted <- count_permutation(real_count_mat = AGP_count_mat, nSample = nSample,
                                  seed=ID)

covariate <- rep(0, nSample)
set.seed(ID)
covariate[sample.int(nSample, nSample*0.5)] <- 1

sample_metadata <- data.frame(X=covariate, 
                              row.names=sprintf("Indv_%d", seq(1, nSample)))

taxa_table <- data.frame(logfold = rep(0, nTaxa),
                         row.names=sprintf("Taxon_%d", seq(1, nTaxa)))

rownames(AGP_permuted) <- sprintf("Taxon_%d", seq(1, nTaxa))
colnames(AGP_permuted) <- sprintf("Indv_%d", seq(1, nSample))

output <- list(count_mat = AGP_permuted,
               sample_metadata = sample_metadata,
               taxa_info =taxa_table)


saveRDS(output, 
        file.path(folder, "data", "null", 
                  sprintf("AGP_simulation_null_%d.rds", ID)))


