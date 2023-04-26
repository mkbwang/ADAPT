
rdirichlet <- function(dir_gamma, nSample, seed=2023){
  set.seed(seed)
  output <- matrix(0, nrow=length(dir_gamma), ncol=nSample)
  for (j in 1:nSample){
    gdist_vec <- rgamma(length(dir_gamma), shape=dir_gamma)
    output[, j] <- gdist_vec / sum(gdist_vec)
  }
  return(output)
}


SimulateCount <- function (
  count_template, nSample = 100,
  diff_prop = 0.1, covariate_type = c("binary", "continuous"),
  grp_ratio = 1, zinf_prop=0.1,
  covariate_eff_mean = 1, covariate_eff_sd = 0, depth_mu = 10000, depth_theta = 5,
  seed=2023
){
  set.seed(seed)
  covariate_type <- match.arg(covariate_type)
  nOTU <- nrow(template_data)

  # permute the counts for each taxon among individuals
  for (j in 1:nOTU){
    template_data[j, ] <- permute(template_data[j, ])
  }
  # estimate dirichlet distribution
  estimated_dirmult <- dirmult::dirmult(t(template_data))
  dir_gamma <- estimated_dirmult$gamma

  # generate relative abundances from dirichlet distribution
  relative_abundance <- rdirichlet(dir_gamma=dir_gamma, nSample=nSample, seed=seed)

  if (covariate_type == 'continuous') {
    X <- rnorm(nSample)
  }
  if (covariate_type == "binary") {
    X <- rep(0, nSample)
    X[sample(1:nSample, round(nSample * grp_ratio / (1+grp_ratio)))] <- 1
  }

  logfoldchange_effect <- rep(0, nOTU) # effect size of log fold change for each taxon
  change_idx <- sample(1:nOTU, round(nOTU*diff_prop)) # taxa which are DA
  # some taxa scale up, others scale down
  change1_idx <- sample(change_idx, round(nOTU*diff_prop*0.5))
  change2_idx <- setdiff(change_idx, change1_idx)
  logfoldchange_effect[change1_idx] <- rnorm(length(change1_idx),
                                          mean=covariate_eff_mean,
                                          sd=covariate_eff_sd)
  logfoldchange_effect[change2_idx] <- rnorm(length(change2_idx),
                                          mean=-covariate_eff_mean,
                                          sd=covariate_eff_sd)


  eta_diff <- cbind(logfoldchange_effect) %*% X
  eta_exp <- exp(eta_diff)
  scaled_relative_abundance <- eta_exp * relative_abundance

  final_relative_abundance <- apply(scaled_relative_abundance, 2,
                                    function(colvec) colvec/sum(colvec))

  nSeq <- rnegbin(nSample, mu = depth_mu, theta = depth_theta)

  otu_tab_sim <- sapply(1:ncol(final_relative_abundance),
                        function(i) rmultinom(1, nSeq[i], final_relative_abundance[,i]))

  # add zero inflation
  nonzero_locs <- which(otu_tab_sim != 0)
  zero_inflation_ids <- sample(nonzero_locs, round(length(nonzero_locs) * zinf_prop),
                               prob=1/otu_tab_sim[nonzero_locs])
  otu_tab_sim[zero_inflation_ids] <- 0

  rownames(otu_tab_sim) <- sprintf("OTU_%d", seq(1, nOTU))
  colnames(otu_tab_sim) <- sprintf("Sample_%d", seq(1, nSample))

  covariate_df <- data.frame(X = X)
  rownames(covariate_df) <- sprintf("Sample_%d", seq(1, nSample))
  foldchange_df <- data.frame(logfold=logfoldchange_effect)
  rownames(foldchange_df) <- sprintf("OTU_%d", seq(1, nOTU))

  simulated_data <- list(otu_tab_sim = otu_tab_sim, metadata = covariate_df,
                         taxa = foldchange_df)

  return(simulated_data)
}



