library(dirmult)
library(phyloseq)
library(dplyr)
library(doParallel)
library(parallelly)
library(fitdistrplus)
library(MASS)

rdirichlet <- function(dir_gamma, nSample, seed=2023){
  set.seed(seed)
  # first generate gamma distribution
  simulation <- sapply(dir_gamma, function(gamma_shape) rgamma(nSample, shape=gamma_shape, rate=1))
  # then normalize to composition
  composition <- apply(simulation, 1, function(sample_count) sample_count/sum(sample_count))
  return(composition)
}

count_permutation <- function(real_count_mat, nSample=100, seed=1){
  # rows are taxa, columns are samples
  if(nSample > ncol(real_count_mat)){
    nSample <- ncol(real_count_mat)
  }
  set.seed(seed)
  permuted_count <- apply(real_count_mat, 1, 
                          function(taxa_counts) taxa_counts[sample.int(ncol(real_count_mat), nSample)]) |> t()
  return(permuted_count)
}


estimate_param <- function(real_count_mat, seed=1){
  # estimate dirichlet multinomial distribution parameter and library sizes based on real data
  # taxa are rows
  set.seed(seed)
  Nsample <- ncol(real_count_mat)
  Ntaxa <- nrow(real_count_mat)
  permuted_count <- count_permutation(real_count_mat = real_count_mat,
                                      nSample = ncol(real_count_mat),
                                      seed=seed)
  
  # negative binomial distribution for sequencing depths
  permuted_depths <- colSums(permuted_count)
  NB_depth_fit <- fitdist(permuted_depths, "nbinom")
  NB_depth_mean <- NB_depth_fit$estimate[2]
  NB_depth_theta <- NB_depth_fit$estimate[1]
  
  # estimate dirichlet multinomial distribution
  dirichlet_param <- dirmult(t(permuted_count))
  # sort the dirichlet parameters
  dirichlet_gamma <- sort(dirichlet_param$gamma, decreasing=TRUE)
  
  result <- list(dirmult_composition=dirichlet_gamma,
                 Seqdepth_mean=NB_depth_mean,
                 Seqdepth_theta=NB_depth_theta)
}



SimulateCount <- function (
  dirmult_composition, seqdepth_mean, seqdepth_theta, nSample = 100,
  grp_ratio = 1, DA_avg_logfold=2, DA_sd_logfold=0.1, DA_proportion=0.2, seed=1, zinf=0.5,
  DA_direction=c("balanced", "unbalanced"), DA_mode=c("abundant", "rare", "mixed") 
){
  DA_direction <- match.arg(DA_direction)
  DA_mode <- match.arg(DA_mode)
  if(DA_proportion >= 0.5){
    stop("The number of DA taxa should not exceed half of all the taxa!")
  }
  # simulate count based on given baseline composition parameter and sequencing depth parameters
  nTaxa <- length(dirmult_composition)
  grp_proportion <- grp_ratio / (grp_ratio + 1) # the sample proportion in two constrasting groups
  
  # sets up the binary covariate
  covariate <- rep(0, nSample)
  set.seed(seed)
  covariate[sample.int(nSample, nSample*grp_proportion)] <- 1

  
  num_DAtaxa <- round(nTaxa * DA_proportion)
  # are all the DA taxa abundant/rare or a mix of both?
  # Note that the dirichlet distribution parameter has been sorted
  set.seed(2023) # fix the taxa which are DA
  if (DA_mode == "mixed"){
    DA_taxa <- sample.int(nTaxa, size=num_DAtaxa)
  } else if(DA_mode == "abundant"){
    DA_taxa <- sample(seq(1, round(nTaxa/4)), size=num_DAtaxa)
  } else if(DA_mode == "rare"){
    DA_taxa <- sample(seq(round(nTaxa*3/4), nTaxa), size=num_DAtaxa)
  }
  
  set.seed(seed)
  log_fold_change_DA <- rnorm(n=num_DAtaxa, 
                              mean=DA_avg_logfold, 
                              sd=DA_sd_logfold) 
  if (DA_direction == "balanced"){
    # Are all the DA taxa changing ine one direction?
    direction <- rep(1, num_DAtaxa)
    direction[sample.int(num_DAtaxa, size=num_DAtaxa/2)] <- -1
    log_fold_change_DA <- log_fold_change_DA * direction
  }
  
  log_fold_change <- rep(0, nTaxa)
  log_fold_change[DA_taxa] <- log_fold_change_DA
  fold_change_mat <- exp(outer(log_fold_change, covariate))
  
  sampled_compositions <- rdirichlet(dirmult_composition, nSample=nSample,
                                     seed=seed)
  # multiply the fold changes and renormalize
  simulated_compositions <- sampled_compositions * fold_change_mat
  simulated_compositions <- apply(simulated_compositions, 2, 
                                  function(comp) comp/sum(comp))
  
  sample_metadata <- data.frame(X=covariate, 
                                row.names=sprintf("Indv_%d", seq(1, nSample)))
  taxa_table <- data.frame(logfold = log_fold_change,
                           row.names=sprintf("Taxon_%d", seq(1, nTaxa)))
  simulation_depths <- rnegbin(n=nSample, mu=seqdepth_mean, theta=seqdepth_theta)
  # simulate count matrix
  cores=parallelly::availableCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  simulated_counts <- foreach(j=1:nSample, .combine=cbind) %dopar%{
    rmultinom(1, size=simulation_depths[j], prob=simulated_compositions[, j])
  }
  stopCluster(cl)
  
  # if there are too few zeros, sample entries and truncate them to zero
  currentzeros <- mean(simulated_counts == 0)
  if (currentzeros < zinf){
    num_extrazeros <- round(length(simulated_counts) * (zinf - currentzeros) )
    nonzero_entries <- which(simulated_counts != 0)
    selected_entries <- sample(nonzero_entries, size=num_extrazeros,
                                   prob=1 / simulated_counts[nonzero_entries])
    simulated_counts[selected_entries] <- 0
  }
  
  
  rownames(simulated_counts) <- sprintf("Taxon_%d", seq(1, nTaxa))
  colnames(simulated_counts) <- sprintf("Indv_%d", seq(1, nSample))
  
  output <- list(count_mat = simulated_counts,
                 sample_metadata = sample_metadata,
                 taxa_info = taxa_table)
  return(output)
}



