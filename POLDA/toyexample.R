library(MASS)
rm(list=ls())

source("POLDA/firth_overdispersed_logit.R")
source("POLDA/overdisperse_GLM.R")

set.seed(2030)
M <- c(rnegbin(n=10, mu=300, theta=5), rnegbin(n=10, mu=20, theta=5))

pi <- rbeta(n=20, shape1=0.5, shape2=50)

Y <- sapply(seq(1, 20), function(j) rbinom(n=1, M[j], pi[j]))
# Y[12] <- 1

X1 <- c(rep(0, 10), rep(1, 10))
X2 <- rep(0, 20)
X2[sample.int(20, 10)] <- 1


X_mat <- cbind(rep(1, 20), X1, X2)
# W <- diag(c(rep(4, 10), rep(9, 10)))

# XWX <- t(X_mat) %*% W %*% X_mat 
# XWX_inv <- chol2inv(chol(XWX))

simulated_df <- data.frame(M=M, Y=Y, X1=X1, X2=X2)

library(dplyr)
simulated_df <- simulated_df %>% arrange(X2)

simulated_df$M_Y <- simulated_df$M - simulated_df$Y

vanilla_LR1 <- glm(cbind(Y, M-Y) ~ X1, 
                  family = binomial("logit"), data = simulated_df)
overdispersed_LR1 <- glm.binomial.disp(vanilla_LR, verbose=FALSE)
overdispersed_LR1$coefficients[2]
sqrt(vcov(overdispersed_LR1)[2,2])

vanilla_LR2 <- glm(cbind(Y, M-Y) ~ X2, 
                   family = binomial("logit"), data = simulated_df)
overdispersed_LR2 <- glm.binomial.disp(vanilla_LR2, verbose=FALSE)
overdispersed_LR2$coefficients[2]
sqrt(vcov(overdispersed_LR2)[2,2])


brglm_fit1 <- brglm(cbind(Y, M-Y) ~ X1, 
                   family = binomial("logit"), data = simulated_df)
brglm_fit1$coefficients[2]
sqrt(vcov(brglm_fit1)[2,2])


brglm_fit2 <- brglm(cbind(Y, M-Y) ~ X2, 
                    family = binomial("logit"), data = simulated_df)
brglm_fit2$coefficients[2]
sqrt(vcov(brglm_fit2)[2,2])


# IRWLS_result <- vanilla_IRWLS(y_vec=Y, m_vec=M, X_mat=X_mat, phi=0.02)
# overdispersed_IRWLS_result <- overdispersed_vanilla_IRWLS(y_vec=Y, m_vec=M, X_mat=X_mat)

overdispersed_firth_IRWLS_result1 <- overdispersed_firth_IRWLS(y_vec=Y, m_vec=M, X_mat=X_mat[, c(1,2)])
beta_hat1 <- overdispersed_firth_IRWLS_result1$beta_hat[2]
sd_beta_hat1 <- sqrt(overdispersed_firth_IRWLS_result1$vcov_beta[2,2])

overdispersed_firth_IRWLS_result2 <- overdispersed_firth_IRWLS(y_vec=Y, m_vec=M, X_mat=X_mat[, c(1,3)])
beta_hat2 <- overdispersed_firth_IRWLS_result2$beta_hat[2]
sd_beta_hat2 <- sqrt(overdispersed_firth_IRWLS_result2$vcov_beta[2,2])


