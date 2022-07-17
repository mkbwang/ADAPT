library(lme4)
library(dplyr)

set.seed(2022)
sythetic_effsizes <- c(0.1, 0.8, 0.3)
covariate_mat <- cbind(1, runif(100, min=-1.5, max=1.5), runif(100, min=-2, max=2))
randeff <- rnorm(100, 3)
logodds <- covariate_mat %*% sythetic_effsizes + randeff
prob <- 1 / (1 + exp(-logodds))
population <- floor(runif(100, min=50, max=100))


# generate long format

ID <- seq(1, 100)
long_ID <- rep(ID, times=population)
covariate_long_mat <- covariate_mat[long_ID, -1]
prob_long <-  rep(prob, times=population)
observed <- rep(0, length(prob_long))
for (j in 1:length(observed)){
  observed[j] <- rbinom(1, 1, prob_long[j])
}
df <- cbind(long_ID, covariate_long_mat, prob_long, observed) %>% as.data.frame()
colnames(df) <- c('ID', 'X1', 'X2', 'Real_Prob', 'Y')


# glmer output
lme4model <- glmer(Y ~ X1 + X2 + (1|ID), data=df, family="binomial",
                   nAGQ=10)


### try with IRWLS

glm_model <- glm(Y ~ X1 + X2, data=df, family="binomial")

beta <- glm_model$coefficients
b <- rep(0, nrow(df))
eta <- as.vector(beta[1] + as.matrix(df[, c('X1', 'X2')]) %*% beta[c(2,3)]) + b
pi <- 1/(1+exp(-eta))
weights <- pi*(1-pi)
df$weight <- weights
u_vec <- eta + (df$Y - pi) / weights
df$U <- u_vec

lmix_model <- lmer(U ~ X1+X2 + (1|ID), data=df, weights = weight, REML=FALSE)
rand_intercepts <- ranef(lmix_model)
## figure out how to iterate the 10 lines above

# some utility functions
# invXRX <- function(xmat, invR){
#   XRX <- t(xmat) %*% invR %*% xmat
#   return(solve(XRX))
# }
#
# proflik <- function(B_scalar, population, pi_vec, u_vec, X_mat){
#   V_vec <- 1 / ( population * pi_vec * (1-pi_vec))
#   R <- B_scalar + V_vec
#   K <- length(V_vec)
#   invmatR <- diag(1/R)
#   XRX_inv <- invXRX(X_mat, invmatR)
#   URX <- t(u_vec) %*% invmatR %*% X_mat
#   URU <- t(u_vec) %*% invmatR %*% u_vec
#   result1 <- as.vector(sum(log(R)))
#   result2 <- as.vector(K * log(URU - URX %*% XRX_inv %*% t(URX)))
#   return(result1 + result2)
# }
#
# ## first IRWLS with regular GLM, this works
# # eta_vec <- observed / population
# # beta_vec <- as.vector(invXRX(covariate_mat, diag(length(eta_vec))) %*% t(covariate_mat) %*% eta_vec)
# # t <- 0
# # repeat{
# #   {
# #     eta_vec <- as.vector(covariate_mat %*% beta_vec)
# #     pi_vec <- 1/(1+exp(-eta_vec))
# #     V_vec <- 1/(population * pi_vec * (1-pi_vec))
# #     u_vec <- eta_vec + V_vec * (observed - population * pi_vec)
# #     new_beta_vec <- as.vector(invXRX(covariate_mat, diag(1/V_vec)) %*% t(covariate_mat) %*%
# #       diag(1/V_vec) %*% u_vec)
# #   };
# #   if(norm(new_beta_vec - beta_vec, type='2') < 0.02){
# #     beta_vec <- new_beta_vec
# #     break
# #   } else{
# #     beta_vec <- new_beta_vec
# #   }
# # }
#
# # confirm that the code above is correct by running the original GLM function
# glm_model <- glm(cbind(Observed, unobserved) ~ u1 + u2, data=testdata, family="binomial")
# sample_num <- nrow(testdata)
# beta_vec <- glm_model$coefficients %>% unname()
# t <- 0 # iteration
# b_vec <- rep(0, 100) # random effects
# eta_vec <- as.vector(covariate_mat %*% beta_vec + b_vec)
# pi_vec <- 1 / (1 + exp(-eta_vec))
# V_vec <- 1 / (population * pi_vec * (1 - pi_vec))
# u_vec <- eta_vec + V_vec * (observed - population * pi_vec)
#
# profliks <- rep(0, 100)
#
# for (j in 1:100){
#    profliks[j]<- proflik(1000+j*100, population, pi_vec, u_vec, covariate_mat)
# }
#
# plot(profliks)
# repeat{
#   {
#     eta_vec <- as.vector(covariate_mat %*% beta_vec + b_vec)
#     pi_vec <- 1 / (1 + exp(-eta_vec))
#     V_vec <- 1 / (population * pi_vec * (1 - pi_vec))
#     u_vec <- eta_vec + V_vec * (observed - population * pi_vec)
#     B_optimize <- nloptr(x0=2,
#                           eval_f=proflik,
#                           lb = 0,
#                           ub = 20,
#                           opts = list("algorithm"="NLOPT_LN_SBPLX",
#                                       "xtol_rel"=1.0e-6),
#                           population = population,
#                           u_vec = u_vec,
#                           pi_vec = pi_vec,
#                           X_mat = covariate_mat) # minimize function
#     B_scalar <- B_optimize$solution
#     R_vec <- B_scalar + V_vec / population
#     invRmat <- diag(1/R_vec)
#     XRX_inv <- invXRX(covariate_mat, invRmat)
#     XRU <- t(covariate_mat) %*% invRmat %*% u_vec
#     beta_new_vec <- as.vector(XRX_inv %*% XRU)
#     URU <- t(u_vec) %*% invRmat %*% u_vec
#     sigma2 <- 1/sample_num * as.vector(URU - t(XRU) %*% XRX_inv %*% XRU)
#     b_vec <- B_scalar/population/R_vec *(u_vec - as.vector(covariate_mat %*% beta_new_vec))
#     t <- t+1
#   };
#   if (t>10){
#     break
#   } else{
#     beta_vec <- beta_new_vec
#   }
# }
#
