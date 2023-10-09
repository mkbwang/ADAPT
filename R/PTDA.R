

#' Pooling Tobit Models for microbiome differential abundance analysis
#' @useDynLib PTDA
#' @param otu_table microbiome abundance table generated from 16S rRNA sequencing or shotgun metagenomic sequencing. All entries are integers
#' @param metadata sample metadata dataframe
#' @param covar the name of the covariate of interest
#' @param adjust the names of confounders to adjust for, NULL if no adjustment
#' @param prevalence_cutoff features whose prevalence are smaller than the cutoff will be excluded from analysis
#' @param depth_cutoff a sample would be discarded if its sequencing depth is smaller than the threshold
#' @param genes_are_rows whether the microbiome genes are on the column or the row of the abundance table
#' @param boot whether to use bootstrap to estimate scaling factors for the test statistics, default TRUE
#' @param boot_replicate number of bootstrap replicates for each gene to estimate scaling factors
#' @param n_boot_gene number of genes to apply bootstrap
#' @param alpha the cutoff of the adjusted p values
#' @importFrom stats optim
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @returns the reference feature set, the identified DA features and the p values for all the features
#' @export
ptda <- function(otu_table, metadata, covar, adjust=NULL,
                 prevalence_cutoff=0.1, depth_cutoff=1000, genes_are_rows=TRUE,
                 boot=TRUE, boot_replicate=1000, n_boot_gene=100, alpha=0.05){

  if(genes_are_rows) otu_table <- t(otu_table)

  # filter out rare taxa and samples with low sequencing depths
  prevalences <- colMeans(otu_table != 0)
  depths <- rowSums(otu_table)
  otu_table_filtered <- otu_table[depths > depth_cutoff, prevalences > prevalence_cutoff]
  gene_names <- colnames(otu_table_filtered)
  refgenes <- gene_names # initially all the taxa are reference taxa(relative abundance)
  # cat("Selecting Reference Set...")
  while(1){
    relabd_output <- count_ratio(otu_table = otu_table_filtered, metadata = metadata,
                                  covar=covar, adjust=adjust, refgenes = refgenes, test_all=FALSE,
                                 boot=F, boot_replicate=boot_replicate, n_boot_gene=n_boot_gene)
    relabd_result <- relabd_output$CR_result
    estimated_effect <- relabd_result$effect
    pvals <- relabd_result$pval
    names(pvals) <- relabd_result$Gene
    names(estimated_effect) <- relabd_result$Gene
    # check distribution of p values
    bumfit <- BUM_fit(pvals)
    loglik <- BUM_llk(bumfit$estim_params, pvals[!is.na(pvals)])
    if (loglik > 2){ # need to continue shrinking reference taxa set
      distance2med <- abs(estimated_effect - median(estimated_effect, na.rm=T))
      sorted_distance <- sort(distance2med)
      ordered_genenames <- names(sorted_distance)
      refgenes <- ordered_genenames[1:(length(ordered_genenames)/2)]
      # cat("Shrink Reference Set...")
    } else{
      break
    }
  }
  # cat("Reference Set Selected, model count ratio of all the genes against reference set ...")
    # after finding the reference set, model the count ratios of all the genes against the reference set
  final_output <- count_ratio(otu_table = otu_table_filtered, metadata = metadata,
                              covar=covar, adjust=adjust, refgenes = refgenes, test_all=T,
                              boot=boot, boot_replicate=boot_replicate, n_boot_gene=n_boot_gene)
  all_CR_results <- final_output$CR_result

  all_pvals <- all_CR_results$pval
  names(all_pvals) <- all_CR_results$Gene
  # find p value cutoff for 0.05 FDR
  all_adjusted_pvals <- p.adjust(all_pvals, method="BH")
  all_CR_results$adjusted_pval <- all_adjusted_pvals

  significant_pvals <- all_adjusted_pvals[all_adjusted_pvals < alpha & !is.na(all_adjusted_pvals)]
  if (length(significant_pvals) > 0){
    DiffGene <- names(significant_pvals)
  } else{
    DiffGene <- c()
  }

  result <- list(Reference_Gene = refgenes,
                 DA_Gene = DiffGene,
                 P_Value = all_CR_results)

  return(result)
}

