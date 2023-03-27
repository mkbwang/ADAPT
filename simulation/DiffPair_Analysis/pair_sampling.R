
rm(list=ls())
## given that the number of pairs could reach millions when the number of taxa reaches thousands,
## It's worth designing a sampling function that samples a subset of all the pairs
## the subset needs to designed such that each taxon is sampled a fixed number of times

pair_subset <- function(num_taxa = 1000, appearance=4){
  stopifnot(appearance < num_taxa - 1)
  results <- matrix(0, ncol=2)
  largest_gap <- num_taxa - 1
  if((appearance %% 2) == 1){
    appearance <- appearance - 1
  }
  gap_upperbound <- appearance / 2
  for (j in 1:gap_upperbound){
    pairs1 <- cbind(seq(1,(num_taxa - j)),
                    seq((1+j),num_taxa))
    pairs2 <- cbind(seq(1, 1+j-1),
                    seq(num_taxa - j+ 1, num_taxa))
    results <- rbind(results, pairs1, pairs2)
  }
  return(results[-1, ])
}

# check counts, confirm that each
final_counts <- rep(0, 100)
trial <- pair_subset(num_taxa = 100, appearance = 6)
for (j in 1:100){
  final_counts[j] <- sum(trial[,1]==j) + sum(trial[,2] == j)
}
stopifnot(all(final_counts == 6))


