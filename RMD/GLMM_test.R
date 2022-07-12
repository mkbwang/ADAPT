library(lme4)
library(dplyr)

set.seed(2022)
sythetic_effsizes <- c(0.1, 0.8, 0.3)
covariate_mat <- cbind(1, runif(100, min=-1.5, max=1.5), runif(100, min=-2, max=2))
randeff <- rnorm(100, 2)
logodds <- covariate_mat %*% sythetic_effsizes + randeff
prob <- 1 / (1 + exp(-logodds))
population <- floor(runif(100, min=50, max=100))
observed <- rep(0, 100)
for (j in 1:100){
  observed[j] <- rbinom(1, population[j], prob[j])
}



testdata <- cbind(seq(1,100), covariate_mat, population, observed) %>% as.data.frame()
colnames(testdata) <- c('ID', 'intercept', 'u1', 'u2', 'Population', 'Observed')
testdata$unobserved <- testdata$Population - testdata$Observed

lme4model <- glmer(cbind(Observed, unobserved) ~ u1 + u2 + (1|ID), data=testdata, family="binomial",
                   nAGQ=10)


### try with IRWLS

# some utility functions
invXRX <- function(xmat, invR){
  XRX <- t(xmat) %*% invR %*% xmat
  eigen_result <- eigen(XRX)
  return (eigen_result$vectors %*% diag(1/eigen_result$values) %*% t(eigen_result$vectors))
}

proflik <- function(B_scalar, population, pi_vec, u_vec, X_mat){
  V_vec <- 1 / ( population * pi_vec * (1-pi_vec))
  R <- B_scalar/population + V_vec
  K <- length(V_vec)
  invmatR <- diag(1/R)
  XRX_inv <- invXRX(X_mat, invmatR)
  URX <- t(u_vec) %*% invmatR %*% X_mat
  URU <- t(u_vec) %*% invmatR %*% u_vec
  result1 <- as.vector(sum(log(R)))
  result2 <- as.vector(K * log(URU - URX %*% XRX_inv %*% t(URX)))
  return(result1 + result2)
}

## first IRWLS with regular GLM, this works
# eta_vec <- observed / population
# beta_vec <- as.vector(invXRX(covariate_mat, diag(length(eta_vec))) %*% t(covariate_mat) %*% eta_vec)
# t <- 0
# repeat{
#   {
#     eta_vec <- as.vector(covariate_mat %*% beta_vec)
#     pi_vec <- 1/(1+exp(-eta_vec))
#     V_vec <- 1/(population * pi_vec * (1-pi_vec))
#     u_vec <- eta_vec + V_vec * (observed - population * pi_vec)
#     new_beta_vec <- as.vector(invXRX(covariate_mat, diag(1/V_vec)) %*% t(covariate_mat) %*%
#       diag(1/V_vec) %*% u_vec)
#   };
#   if(norm(new_beta_vec - beta_vec, type='2') < 0.02){
#     beta_vec <- new_beta_vec
#     break
#   } else{
#     beta_vec <- new_beta_vec
#   }
# }

# confirm that the code above is correct by running the original GLM function
glm_model <- glm(cbind(Observed, unobserved) ~ u1 + u2, data=testdata, family="binomial")
sample_num <- nrow(testdata)
beta_vec <- glm_model$coefficients %>% unname()
t <- 0 # iteration
b_vec <- rep(0, 100) # random effects
eta_vec <- as.vector(covariate_mat %*% beta_vec + b_vec)
pi_vec <- 1 / (1 + exp(-eta_vec))
V_vec <- 1 / (population * pi_vec * (1 - pi_vec))
u_vec <- eta_vec + V_vec * (observed - population * pi_vec)

results_1 <- rep(0, 100)
results_2 <- rep(0, 100)

for (j in 1:100){
  currentresult <- proflik(j*0.001, population, pi_vec, u_vec, covariate_mat)
  results_1[j] <- currentresult[[1]]
  results_2[j] <- currentresult[[2]]
}

repeat{
  {
    eta_vec <- as.vector(covariate_mat %*% beta_vec + b_vec)
    pi_vec <- 1 / (1 + exp(-eta_vec))
    V_vec <- 1 / (population * pi_vec * (1 - pi_vec))
    u_vec <- eta_vec + V_vec * (observed - population * pi_vec)
    B_optimize <- nloptr(x0=2,
                          eval_f=proflik,
                          lb = 0,
                          ub = 20,
                          opts = list("algorithm"="NLOPT_LN_SBPLX",
                                      "xtol_rel"=1.0e-6),
                          population = population,
                          u_vec = u_vec,
                          pi_vec = pi_vec,
                          X_mat = covariate_mat) # minimize function
    B_scalar <- B_optimize$solution
    R_vec <- B_scalar / population + V_vec
    invRmat <- diag(1/R_vec)
    XRX_inv <- invXRX(covariate_mat, invRmat)
    XRU <- t(covariate_mat) %*% invRmat %*% u_vec
    beta_new_vec <- as.vector(XRX_inv %*% XRU)
    URU <- t(u_vec) %*% invRmat %*% u_vec
    sigma2 <- 1/sample_num * as.vector(URU - t(XRU) %*% XRX_inv %*% XRU)
    b_vec <- B_scalar/population/R_vec *(u_vec - as.vector(covariate_mat %*% beta_new_vec))
    t <- t+1
  };
  if (t>10){
    break
  } else{
    beta_vec <- beta_new_vec
  }
}

