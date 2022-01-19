library(testthat)
library(DiffRatio)
library(extraDistr)

covar1 = rep(c(0, 1), each=21)
covar2 = rep(seq(-10, 10), 2)

# taxa mean
tax1_mean = round(exp(8.5+covar1 + 0.15*covar2))
tax2_mean = round(exp(9.5-covar1 + 0.2*covar2))
tax3_mean = round(rep(exp(8.5), 42))
tax4_mean = round(rep(exp(9), 42))
tax5_mean = round(rep(exp(10), 42))

allconditions = cbind(covar1, covar2, tax1_mean, tax2_mean, tax3_mean, tax4_mean, tax5_mean)
rep_conditions <- allconditions[rep(seq_len(nrow(allconditions)), each = 20), ]
population = matrix(0, 840, 5) # total abundance generated using negative binomial distribution
counts = matrix(0, 840, 5) # sampled from multihypergeometric distribution
depth <- round(rnorm(840, mean=15000, sd=1000)) # sequencing depth

for (id in 1:nrow(population)){
  newsample <- rep(0, 5) # a new sample population
  for (pos in 1:5){
    newsample[pos] <- rnbinom(n=1, mu=rep_conditions[id, pos+2], size=1e5)
  }
  population[id,] <- newsample
  counts[id, ]<- rmvhyper(1, newsample, depth[id])
}

colnames(population) <- c('Tax1_pop', 'Tax2_pop', 'Tax3_pop', 'Tax4_pop', 'Tax5_pop')
colnames(counts) <- c('Tax1_count', 'Tax2_count', 'Tax3_count', 'Tax4_count', 'Tax5_count')

testdata <- cbind(rep_conditions, population, counts, depth)
reg_result <- dfr(data=testdata, covar=c("covar1", "covar2"), tpair=c("Tax1_count", "Tax2_count"))

test_that("Logistic Regression working", {
  expect_equal(reg_result$coefficients[['covar1']], 2, tolerance=1e-2)
  expect_equal(reg_result$coefficients[['covar2']], -0.05, tolerance=1e-2)
})
