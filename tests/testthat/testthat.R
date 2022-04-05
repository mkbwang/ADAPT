library(testthat)
library(DiffRatio)

beta0 <- -0.5 # baseline average log odds
beta1 <- 1.5 # average effect for each individual

# add in random effects for each individual
baseline_logodds <- rnorm(30, mean=beta0, sd=0.1)
effect_logodds <- rnorm(15, mean=beta1, sd=0.1)

# log odds of two taxa of interest for 30 individuals
logodds <- baseline_logodds + c(rep(0, 15), effect_logodds)
proportions <- cbind(exp(logodds)/(1+exp(logodds)), 1/(1+exp(logodds)))

# take into account several other ignorant taxa
full_proportions <- cbind(proportions * 0.4,
                          t(replicate(30, c(0.15, 0.2, 0.25))))

# read depths for each individual
read_depths <- runif(30, min=7000, max=13000)

# counts for each individual
count_mat <- matrix(0, nrow=30, ncol=5)
for (i in 1:nrow(count_mat)){
  count_mat[i,] <- rmultinom(n=1, read_depths[i], full_proportions[i, ])
}


colnames(count_mat) <- c("taxa1", "taxa2", "taxa3", "taxa4", "taxa5")
count_mat <- otu_table(count_mat, taxa_are_rows = FALSE)
X <- c(rep(0, 15), rep(1, 15))
sample_info <- sample_data(data.frame(X))
test_data <- phyloseq(count_mat, sample_info)

reg_result1 <- dfr(data=test_data, covar=c("X"), tpair=c("taxa1", "taxa2"))

test_that("Logistic Regression working", {
  expect_equal(reg_result1$coefficients[[1]], -0.5, tolerance=1e-1)
  expect_equal(reg_result1$coefficients[[2]], 1.5, tolerance=1e-1)
})
