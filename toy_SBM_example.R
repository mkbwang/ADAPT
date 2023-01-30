# SBM toy example, illustrate the essentials of blockmodels package

n <- 12 # nodes
Z <- matrix(0, nrow=n, ncol=2)
Z[1:6, 1] <- 1; Z[7:12, 2] <- 1
# Z<-diag(Q)%x%matrix(1,npc,1)
P<-matrix(c(0.6, 0.9, 0.9, 0.2), nrow=2)
P[lower.tri(P)]<-t(P)[lower.tri(P)]
set.seed(20)
M<-1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z)) ## adjacency matrix
M[lower.tri(M)]<-t(M)[lower.tri(M)]
diag(M) <- 0


# single block probability

singleblock_prob <- sum(M)/(nrow(M)^2 - nrow(M))

## calculate eigen vectors and carry out spectral clustering

error <- M - singleblock_prob
W<- error %*% t(error)
W<-1/(1+exp(-W/sd(W)))
D<- diag(1/sqrt(rowSums(W)))
L<- D %*% W %*% D

eig_decomp <- eigen(L, symmetric = TRUE)

## kmeans to get initial categorical membership

eigvecs <- eig_decomp$vectors[, c(1,2)]
centroids <- matrix(0, 2, 2)
dists <- as.matrix(dist(eigvecs))
choosen <- which(dists==max(dists),arr.ind = TRUE)[1,c('row','col')]
centroids <- eigvecs[choosen,]

starting_point <- kmeans(eigvecs, centers=centroids, iter.max=20)

init_membership <- matrix(0, nrow=nrow(eigvecs), ncol=2)
init_membership[starting_point$cluster == 1, 1] <- 1
init_membership[starting_point$cluster == 2, 2] <- 1


# cap the membership probability value so that it is always nonzero
cap_membership <- function(memmat, upper_bound=1-1e-6, lower_bound=1e-6){
  memmat[memmat > upper_bound] <- upper_bound
  memmat[memmat < lower_bound] <- lower_bound
  return(memmat/rowSums(memmat))
}


# initialize all the quantities for EM
membership <- cbind(colMeans(M) , 1-colMeans(M))
membership <- cap_membership(init_membership) # membership probability matrix
adj_ZD <- M # adjacency matrix (diagonal zero)


alpha_update <- function(membermat){
  colSums(membermat) / sum(membermat)
}

entropy_update <- function(membermat){
  -sum(membermat * log(membermat))
} # calculate sum of entropy


pi_update <- function(adjmat, membermat){
  onesZD <- matrix(1, nrow=nrow(membermat), ncol=nrow(membermat)) -
    diag(nrow(membermat))
  return((t(membermat) %*% adjmat %*% membermat) /
           (t(membermat) %*% onesZD %*% membermat))
} # maximization step, update block probability

pseudolik <- function(adjmat, membermat, alphavec, pimat){
  prior <- sum(membermat %*% log(alphavec))
  onesZD <- matrix(1, nrow=nrow(membermat), ncol=nrow(membermat)) -
    diag(nrow(membermat))
  posterior <- 0.5 * sum((log(pimat) - log(1-pimat)) * (t(membermat) %*% adjmat %*% membermat)) +
    0.5 * sum(log(1-pimat) * (t(membermat) %*% onesZD %*% membermat))

  return(prior + posterior)
}


membership_update <- function(adjmat, membermat, alpha_vec, pimat){
  one_minus_adj_ZD <- 1 - adjmat - diag(nrow(membermat))
  repeat{
    lZ <- t(replicate(nrow(membermat), log(alpha_vec))) # prior
    lZ <- lZ + adjmat %*% membermat %*% log(pimat) +
      one_minus_adj_ZD %*% membermat %*% log(1-pimat) # plus posterior
    lZ <- lZ - rowMaxs(lZ) # shift
    new_Z <- exp(lZ) / rowSums(exp(lZ)) # normalize
    new_membermat <- cap_membership(new_Z) # cap extreme high/low probability
    step_size <- max(abs(membermat - new_membermat))
    membermat <- new_membermat
    if (step_size < 1e-2){
      break
    }
  }
  return(membermat)
} # E step, update


alpha_vec <- colSums(membership) / sum(membership) # overall group proportion
pi_mat <- pi_update(adj_ZD, membership) # block probability
entropy <- entropy_update(membership)
old_ELBO <- entropy + pseudolik(adj_ZD, membership, alpha_vec, pi_mat)
ELBO <- old_ELBO

# EM iteration
repeat{
  ## E step
  membership <- membership_update(adj_ZD, membership, alpha_vec, pi_mat)
  ## M steps
  alpha_vec <- colSums(membership) / sum(membership) # overall group proportion
  pi_mat <- pi_update(adj_ZD, membership) # block probability

  entropy <- entropy_update(membership)
  ELBO <- entropy + pseudolik(adj_ZD, membership, alpha_vec, pi_mat)
  print(ELBO)
  stepsize <- abs(ELBO - old_ELBO)
  old_ELBO <- ELBO
  if (stepsize < 1e-3){
    break
  }
}

