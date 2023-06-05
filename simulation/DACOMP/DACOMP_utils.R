library(pROC)
library(dacomp)

# dacomp pipeline
pipeline <- function(AGP_data, alpha=0.05){
  sample_metadata <- AGP_data$sample_metadata
  count_mat <- AGP_data$count_mat
  
  # find a reference set
  selected_references <- dacomp.select_references(
    X = t(count_mat), verbose = F,
    minimal_TA=10)
  
  updated_selected_references <- dacomp.validate_references(X = t(count_mat), #Counts
                                                            Y = sample_metadata$X, #Traits
                                                            ref_obj = selected_references, #reference checked, must be object from dacomp.select_references(...)
                                                            test = DACOMP.TEST.NAME.WILCOXON, #Test used for checking, can be same test used in dacomp.test(...)
                                                            Q_validation = 0.1, #FDR level for checking
                                                            Minimal_Counts_in_ref_threshold = 1, #reference taxa will must include at least this number of reads
                                                            Reduction_Factor = 0.1, #multiplicative factor used for lowering the threshold for the number of reads required in reference taxa at each iteration
                                                            NR_perm = 1000, #number of permutations used for testing. should be at least 1/(Q_validation/ncol(X))
                                                            Verbose = F)
  
  # run test based on reference set
  test_result = dacomp.test(X = t(count_mat), #counts data
                            y = sample_metadata$X, #phenotype in y argument
                            # obtained from dacomp.select_references(...):
                            ind_reference_taxa = updated_selected_references,
                            test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                            disable_DSFDR = T,
                            q=alpha,
                            verbose = F)
  
  rejected_BH <- which(test_result$p.values.test.adjusted <= alpha)
  
  output <- list(raw_pval = test_result$p.values.test,
                 adjusted_pval = test_result$p.values.test.adjusted,
                 DA_taxa = rejected_BH,
                 ref_taxa = updated_selected_references)
  return(output)
  
}


# evaluate performance
evaluation <- function(taxa_truth, dacomp_result, nullcase=FALSE){
  result <- NULL
  if (nullcase){
    raw_pval <- dacomp_result$raw_pval
    FPR <- mean(raw_pval < 0.05, na.rm=TRUE)
    result <- list(FPR=FPR)
  } else{
    raw_pval <- dacomp_result$raw_pval
    adjusted_pval <- dacomp_result$adjusted_pval
    all_taxa <- rownames(taxa_truth)
    dacomp_DA_taxa <- sprintf("Taxon_%d", dacomp_result$DA_taxa)
    dacomp_reftaxa <- sprintf("Taxon_%d", dacomp_result$ref_taxa)
    true_DA_taxa <- all_taxa[taxa_truth$logfold != 0]
    DA_truth_binary <- taxa_truth[!is.na(raw_pval), ] == 0
    roc_obj <- roc(DA_truth_binary, raw_pval[!is.na(raw_pval)])
    
    check_reftaxa <- dacomp_reftaxa %in% true_DA_taxa
    reftaxa_error <- mean(check_reftaxa)
    check_DAtaxa <- dacomp_DA_taxa %in% true_DA_taxa
    FDR <- 1 - mean(check_DAtaxa)
    Power <- sum(check_DAtaxa) / length(true_DA_taxa)
    result <- list(reftaxa_error = reftaxa_error,
                   FDR = FDR,
                   Power=Power,
                   AUROC = roc_obj$auc)
  }
  
  return(result)
}


