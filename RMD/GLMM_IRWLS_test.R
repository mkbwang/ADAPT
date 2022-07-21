library(lme4)
library(dplyr)

source('RMD/preprocess.R')
set.seed(2022)
# sythetic_effsizes <- c(1.2, -1.3)
# covariate_mat <- cbind(rbinom(100, 1, 0.5), rbinom(100, 1, 0.5))

### random effects generated from two Gaussian mixtures
# mixcomp1 <- rnorm(100, mean=-2, sd=0.4)
# mixcomp2 <- rnorm(100, mean=2, sd=0.4)
# selector <- rbinom(100, 1, prob=0.5)
# randeff <- rep(0, 100)
# randeff[selector == 0] <- mixcomp1[selector == 0]
# randeff[selector == 1] <- mixcomp2[selector == 1]


# logodds <- 0.05 + covariate_mat %*% sythetic_effsizes + randeff
# prob <- 1 / (1 + exp(-logodds))
# population <- floor(runif(100, min=10, max=25))
#

# generate data

# ID <- seq(1, 100)
# positive <- rep(0, length(ID))
# negative <- rep(0, length(ID))
# for (j in 1:length(ID)){
#   positive[j] <- rbinom(1, population[j], prob[j])
#   negative[j] <- population[j] - positive[j]
# }
# df <- cbind(ID, covariate_mat, positive, negative) %>% as.data.frame()
# colnames(df) <- c('ID', 'X1', 'X2', 'Positive', 'Negative')
# head(df)

# transform data from short to long
short2long <- function(shortdf, count_cols){

  covar_df <- shortdf[, !names(shortdf) %in% count_cols]

  cnames <- c(colnames(covar_df), 'Y')
  longdf <- data.frame(matrix(ncol = length(cnames), nrow=0))
  colnames(longdf) <- cnames

  counts_df <- shortdf[, count_cols]
  freq <- rep(0, nrow(counts_df))
  for (j in 1:nrow(counts_df)){
    pos_count <- counts_df[j, count_cols[1]]
    neg_count <- counts_df[j, count_cols[2]]
    if (pos_count + neg_count == 0) next
    freq[j] <- pos_count + neg_count
    Y <- c(rep(1, pos_count), rep(0, neg_count))
    repeat_covars <- covar_df[rep(j, pos_count+neg_count), ]
    longdf <- rbind(longdf, cbind(repeat_covars, Y))
  }
  freq <- freq[freq!=0]
  return(list(longdf, freq))
}

# long_df <- short2long(df, c('Positive', 'Negative'))
# df$ID <- as.factor(df$ID)
# colnames(long_df) <- c('ID', 'X1', 'X2', 'Y')


# glmer output
# glmermodel <- glmer(cbind(Positive, Negative) ~ X1 + X2 + (1|ID), data=df, family="binomial",
#                    nAGQ = 0)

### try with IRWLS

IRWLS <- function(taxa_pairs, covariates, level,  data){
  # generate long format
  clean_data <- short2long(data, taxa_pairs)
  long_df <- clean_data[[1]]
  freq <- clean_data[[2]]
  long_df[[level]] <- as.factor(long_df[[level]])

  # generate regression formula
  covariates_sum <- paste(covariates, collapse='+')
  reffect <- sprintf('(1|%s)', level)
  regfml_glm <- formula(sprintf("%s %s", 'Y ~', covariates_sum))
  regfml_lmm <- formula(sprintf("%s %s + %s", 'U ~', covariates_sum, reffect))

  # glm
  glm_model <- glm(regfml_glm, data=long_df, family="binomial")

  beta <- glm_model$coefficients
  b <- rep(0, length(unique(long_df[[level]])))
  t <- 0

  # Iteratively weighted LMM
  repeat{
    {
      eta <- beta[1] + rep(b, freq) +
        as.vector(as.matrix(long_df[, covariates]) %*% beta[-1])
      pi <- 1/(1+exp(-eta))
      weights <- pi*(1-pi)
      long_df$weight <- weights
      u_vec <- eta + (long_df$Y - pi) / weights
      long_df$U <- u_vec
      lmix_model <- lmer(regfml_lmm, data=long_df, weights = weight, REML=FALSE)
      new_beta <- fixef(lmix_model) %>% unname()
      b <- ranef(lmix_model)[[1]][[1]]
      t <- t+1
    };
    if(norm(new_beta - beta, type='2') < 0.000001){
      beta <- new_beta
      final_model <- lmix_model
      break;
    } else{
      beta <- new_beta
    }
  }
  return(final_model)
}

# final_model <- IRWLS(c('Positive', 'Negative'), c('X1', 'X2'), 'ID', df)


get_pval <- function(tvar){
  0.5 - abs(0.5 - pnorm(tvar))
}

# beta_glmer <- fixef(glmermodel) %>% unname()
# var_glmer <- vcov(glmermodel) %>% diag()
# t_glmer <- beta_glmer / sqrt(var_glmer)
# pval_glmer <- get_pval(t_glmer)
#
# beta_IRWLS <- fixef(final_model)
# var_IRWLS <- vcov(final_model) %>% diag()
# t_IRWLS <- beta_IRWLS / sqrt(var_IRWLS)
# pval_IRWLS <- get_pval(t_IRWLS)



## try on real data, especially when glmer fails, first find two pairs
pvals_glmer = read.csv('RMD/CAARS_Model_Summary/Pvals.csv')
pvals_mat <- pvals_glmer[,c('GLMM_PIRLS', 'GLMM_LA_NM', 'GLMM_LA_BOBYQA', 'GLMM_GQ_NM',
                            'GLMM_GQ_BOBYQA')] %>% as.matrix()
pval_var <- apply(pvals_mat, 1, var)
tpair_pval_rank <- length(pval_var) + 1 - rank(pval_var)
rowids <- which(tpair_pval_rank %in% seq(1, 10))
selected_pairs <- pvals_glmer[rowids, ] # two pairs of taxa to demonstrate


first_pair <- c(selected_pairs$Taxa1[1], selected_pairs$Taxa2[1])
second_pair <- c(selected_pairs$Taxa1[2], selected_pairs$Taxa2[2])

### download data

folder = 'RMD/CAARS_data'
load(file.path(file.path(folder, 'CAARS_processed_GENUS.Rdata')))
processed_data <- preprocess(CAARS.data.genus, 'SAMPLE_ID', "asthma")


filtered_count <- processed_data$feature_table
filtered_metadata <- processed_data$meta_data
struc_zero <- processed_data$structure_zeros

covariates <- filtered_metadata[, c('SAMPLE_ID', 'asthma')]
df1 <- cbind(covariates, t(filtered_count[first_pair, ]))
glmermodel1 <- glmer(cbind(g__Saccharimonadales, g__Moraxella) ~ asthma + (1|SAMPLE_ID),
                     data=df1, family="binomial")
IRWLSmodel1 <- IRWLS(c('g__Saccharimonadales', 'g__Moraxella'), c('asthma'), 'SAMPLE_ID', df1)

df2 <- cbind(covariates, t(filtered_count[second_pair, ]))
glmermodel2 <- glmer(cbind(g__Streptobacillus, g__Kingella) ~ asthma + (1|SAMPLE_ID),
                     data=df2, family="binomial")
IRWLSmodel2 <- IRWLS(c('g__Streptobacillus', 'g__Kingella'), c('asthma'), 'SAMPLE_ID', df2)

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
