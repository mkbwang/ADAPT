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

count_df <- as.data.frame(count_mat)
colnames(count_df) <- c("taxa1", "taxa2", "taxa3", "taxa4", "taxa5")
count_df$id <- seq(1, 30)
count_df$X <- c(rep(0, 15), rep(1, 15))


# model <- glmer(cbind(taxa1, taxa2) ~ X + (1|id), data=count_df, family="binomial")
# model2 <- glm(cbind(taxa1, taxa2) ~ X , data=count_df, family="binomial")
# testdata <- cbind(rep_conditions, population, counts, depth)
reg_result1 <- dfr(data=count_df, covar=c("X"), tpair=c("taxa1", "taxa2"))
reg_result2 <- dfr(data=count_df, covar=c("X"), tpair=c("taxa1", "taxa2"), indveff="id")

test_that("Logistic Regression working", {
  expect_equal(reg_result1$coefficients[[1]], -0.5, tolerance=1e-1)
  expect_equal(reg_result1$coefficients[[2]], 1.5, tolerance=1e-1)
  expect_equal(reg_result2@beta[1], -0.5, tolerance=1e-1)
  expect_equal(reg_result2@beta[2], 1.5, tolerance=1e-1)
})
