# SBM_utils
library(MatrixGenerics)

## prevent probability estimate from being exactly zero during the iteration
cap_membership <- function(memmat, gap=1e-5){
  upper_bound <- 1 - gap/nrow(memmat)
  lower_bound <- gap/nrow(memmat)
  memmat[memmat > upper_bound] <- upper_bound
  memmat[memmat < lower_bound] <- lower_bound
  return(memmat/rowSums(memmat))
}

## E step update of node(taxon) membership (DA/non-DA)
membership_update <- function(adjmat, membermat, pi_vec, theta_mat){
  one_minus_adj_ZD <- 1 - adjmat - diag(nrow(membermat))
  lZ <- t(replicate(nrow(membermat), log(pi_vec))) # prior
  lZ <- lZ + adjmat %*% membermat %*% log(theta_mat) +
    one_minus_adj_ZD %*% membermat %*% log(1-theta_mat) # plus posterior
  lZ <- lZ - rowMaxs(lZ) # shift
  new_Z <- exp(lZ) / rowSums(exp(lZ)) # normalize
  new_membermat <- cap_membership(new_Z) # cap extreme high/low probability
  return(new_membermat)
}


## M step, update the overall membership allocation (DA/non-DA proportion)
pi_update <- function(membermat){
  colSums(membermat) / sum(membermat)
}

## M step, update the block connection probability
theta_update <- function(adjmat, membermat, restriction = FALSE){

  onesZD <- matrix(1, nrow=nrow(membermat), ncol=nrow(membermat)) -
    diag(nrow(membermat))
  updated_theta_mat <- (t(membermat) %*% adjmat %*% membermat) /
    (t(membermat) %*% onesZD %*% membermat)
  ## we want to see the connection probability in off-diagonal blocks being bigger than diagonal terms
  if (restriction & (updated_theta_mat[1,2] < max(diag(updated_theta_mat)))){
    updated_theta_mat[1,2] <- max(diag(updated_theta_mat))
    updated_theta_mat[2,1] <- max(diag(updated_theta_mat))
  }

  return(updated_theta_mat)
} # maximization step, update block probability


## ELBO calculation
entropy_update <- function(membermat){
  -sum(membermat * log(membermat))
} # entropy


pseudolik <- function(adjmat, membermat, pi_vec, theta_mat){
  prior <- sum(membermat %*% log(pi_vec)) # prior log likelihood
  onesZD <- matrix(1, nrow=nrow(membermat), ncol=nrow(membermat)) -
    diag(nrow(membermat))
  posterior <- 0.5 * sum((log(theta_mat) - log(1-theta_mat)) *
                           (t(membermat) %*% adjmat %*% membermat)) +
    0.5 * sum(log(1-theta_mat) * (t(membermat) %*% onesZD %*% membermat))

  return(prior + posterior)
} # expected log likelihood of the adjacency matrix and the membership(DA) probability



## complete EM function
complete_EM <- function(adj, estim_membership, restriction=FALSE){

  membership_iterlist <- list()
  pi_vec <- pi_update(estim_membership) # overall group proportion

  theta_mat <- theta_update(adj, estim_membership) # block probability
  ## we want to see the connection probability in off-diagonal blocks being bigger than diagonal terms
  if (restriction & (theta_mat[1,2] < max(diag(theta_mat)))){
    theta_mat[1,2] <- max(diag(theta_mat))
    theta_mat[2,1] <- max(diag(theta_mat))
  }

  entropy <- entropy_update(estim_membership)
  old_ELBO <- entropy + pseudolik(adj, estim_membership, pi_vec, theta_mat)
  ELBO <- old_ELBO
  counter <- 0 # iteration counter

  # EM iteration
  repeat{
    if(counter == 0){
      # cat("Initial group membership: \n", round(estim_membership[, 1], digits=4), "\n")
      cat("Initial ELBO: ", ELBO, "\n")
      membership_iterlist <- append(membership_iterlist, estim_membership)
    }

    ## E step
    estim_membership <- membership_update(adj, estim_membership, pi_vec, theta_mat)
    ## M steps
    pi_vec <- pi_update(estim_membership) # overall group proportion
    theta_mat <- theta_update(adj, estim_membership, restriction=restriction) # block probability
    entropy <- entropy_update(estim_membership)
    ELBO <- entropy + pseudolik(adj, estim_membership, pi_vec, theta_mat)
    counter <- counter + 1
    # cat("Group membership after iteration ", counter, ": \n", round(estim_membership[, 1], digits=4), '\n')
    cat("ELBO after iteration ", counter, ": ", ELBO, '\n')

    stepsize <- abs(ELBO - old_ELBO)
    old_ELBO <- ELBO
    if (stepsize < 1e-5){
      break
    }
  }
  return(list(membership=estim_membership, theta=theta_mat, ELBO = ELBO))
}


## membership initialization using spectral clustering
spectral_clustering <- function(adj){

  ## take residual
  residual <- adj - sum(adj)/(nrow(adj)^2 - nrow(adj))
  W<- residual %*% t(residual)
  W<-1/(1+exp(-W/sd(W)))
  D<- diag(1/sqrt(rowSums(W)))
  L<- D %*% W %*% D
  eig_decomp <- eigen(L, symmetric = TRUE)

  ## kmeans to get initial categorical membership

  coordinates <- eig_decomp$vectors[, c(1,2)]
  centroids <- matrix(0, 2, 2)
  dists <- as.matrix(dist(coordinates))
  chosen <- which(dists==max(dists),arr.ind = TRUE)[1,c('row','col')]
  centroids <- coordinates[chosen,]
  starting_point <- kmeans(coordinates, centers=centroids, iter.max=20)
  init_membership <- matrix(0, nrow=nrow(coordinates), ncol=2)
  init_membership[starting_point$cluster == 1, 1] <- 1
  init_membership[starting_point$cluster == 2, 2] <- 1

  return(cap_membership(init_membership))

}


## membership probability initialization using degree numbers of each node(taxon)
degree_init <- function(adj){
  avg_degree <- colMeans(adj) * nrow(adj) / (nrow(adj) - 1)
  unname(cbind(avg_degree, 1 - avg_degree))
}


## flat membership initialization, equivalent to single block
flat_init <- function(adj){
  matrix(0.5, nrow=nrow(adj), ncol=2)
}





