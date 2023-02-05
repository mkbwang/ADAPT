## ----setup, include=FALSE----------------------------------------------------------------
remove(list=ls())

source('simulation/DiffTaxon_ID/SBM_utils.R')
library(blockmodels)

## ----------------------------------------------------------------------------------------
n <- 12 # nodes
real_membership <- matrix(0, nrow=n, ncol=2)
real_membership[1:6, 1] <- 1; real_membership[7:12, 2] <- 1


## easy example (distinguishable block probabilities)
real_theta_easy<-matrix(c(0.6, 0.9, 0.9, 0.2), nrow=2)
set.seed(20)
adj_easy <-1*(matrix(runif(n*n),n,n)<real_membership %*%
                real_theta_easy %*% t(real_membership)) ## adjacency matrix
adj_easy[lower.tri(adj_easy)]<-t(adj_easy)[lower.tri(adj_easy)]
diag(adj_easy) <- 0
adj_easy


## hard example (undistinguishable block probabilities)
real_theta_hard<-matrix(c(0.5, 0.65, 0.65, 0.35), nrow=2)
set.seed(20)
adj_hard <-1*(matrix(runif(n*n),n,n)<real_membership %*%
                real_theta_hard %*% t(real_membership)) ## adjacency matrix
adj_hard[lower.tri(adj_hard)]<-t(adj_hard)[lower.tri(adj_hard)]
diag(adj_hard) <- 0
adj_hard


## ----------------------------------------------------------------------------------------
init_membership_sc <- spectral_clustering(adj_easy)
result_sc <- complete_EM(adj_easy, init_membership_sc)
membership_sc <- result_sc$membership

library(blockmodels)
result_package <- BM_bernoulli(membership_type="SBM_sym", adj=adj_easy,
                                          explore_min=2, explore_max=2, plotting='', verbosity=6)
result_package$estimate()
membership_package_sc <- result_package$memberships[[2]]$Z
difference <- membership_package_sc - membership_sc
max(abs(difference))

## ----------------------------------------------------------------------------------------
init_membership_degree <- cap_membership(degree_init(adj_easy))
result_degree <- complete_EM(adj_easy, init_membership_degree)


## ----------------------------------------------------------------------------------------
init_membership_flat <- cap_membership(flat_init(adj_easy))
result_flat <- complete_EM(adj_easy, init_membership_flat)


## ----------------------------------------------------------------------------------------
cat("Block connection probability after spectral clustering initialization is:\n")
result_sc$theta
cat("Block connection probability after degree-based initialization is:\n")
result_degree$theta
cat("Block connection probability after flat initialization is:\n")
result_flat$theta


## ----------------------------------------------------------------------------------------
init_membership_sc <- spectral_clustering(adj_hard)
result_sc <- complete_EM(adj_hard, init_membership_sc)
result_package <- BM_bernoulli(membership_type="SBM_sym", adj=adj_hard,
                               explore_min=2, explore_max=2, plotting='', verbosity=6)
result_package$estimate()

membership_sc <- result_sc$membership
membership_package_sc <- result_package$memberships[[2]]$Z
difference <- membership_package_sc - membership_sc
max(abs(difference))

## ----------------------------------------------------------------------------------------
init_membership_degree <- cap_membership(degree_init(adj_hard))
result_degree <- complete_EM(adj_hard, init_membership_degree)


## ----------------------------------------------------------------------------------------
init_membership_flat <- cap_membership(flat_init(adj_hard))
result_flat <- complete_EM(adj_hard, init_membership_flat)


## ----------------------------------------------------------------------------------------
cat("Block connection probability after spectral clustering initialization is:\n")
result_sc$theta
cat("Block connection probability after degree-based initialization is:\n")
result_degree$theta
cat("Block connection probability after flat initialization is:\n")
result_flat$theta

