rm(list=ls())
# load data

data_withDA <- readRDS(file.path('simulation', 'data', 'PoiGamma', 'simple_simulation_withDA2.rds'))
truth <- data_withDA$mean.eco.abn
glm_result <- read.csv(file.path('simulation', 'glmdisp_result', 'glmdisp_simple_result_withDA2.csv'))


Z_score <- glm_result$effect / glm_result$SE
Z_score <- c(0, Z_score) # assume that we have observed theta_1 to be zero
num_taxa <- 100
allpairs <- combn(seq(1, num_taxa), m=2)
design_matrix <- matrix(0, nrow=nrow(glm_result), ncol=100)
for (j in 1:ncol(allpairs)){
  design_matrix[j, allpairs[1, j]] <- 1
  design_matrix[j, allpairs[2, j]] <- -1
}
design_matrix <- rbind(c(1, rep(0, num_taxa - 1)), design_matrix)

# assume that all the thetas have the same variance sigma^2 and all the z scores given
# thetas have variance 1/4*sigma^2

V_mat <- design_matrix %*% t(design_matrix)
V_mat <- V_mat + diag(1/4, nrow=4951, ncol=4951)
V_mat[1,1] <- 1

# this inverse takes a long time
# need to improve (e.g. only fit a subset of logistic regressions)
ptm <- proc.time()
V_inv <- chol2inv(chol(V_mat))
proc.time() - ptm


XtVX <- t(design_matrix) %*% V_inv %*% design_matrix
XtVX_inv <- chol2inv(chol(XtVX))

# generalized least square estimate
tau_hat <- XtVX_inv %*% t(design_matrix) %*% V_inv %*% Z_score
# residual
resid <- Z_score - design_matrix %*% tau_hat
# sigma^2 estimate (MLE)
sigma2_hat <- t(resid) %*% V_inv %*% resid / nrow(resid)
covar_tau <- drop(sigma2_hat) * XtVX_inv
stderror_tau <- sqrt(diag(covar_tau))

# visualize the distribution of estimated taus for each taxon
library(ggplot2)
strength <- data.frame(estimate = tau_hat)
ggplot(strength, aes(x=estimate)) +
  geom_histogram(color="black", fill="white", binwidth=1)+
  xlab("Estimated Strength Correction")


# we notice that the mode(peak of density) isn't at zero. This is because the selected theta_1 isn't actually zero.
# No matter, we can shift all the tau estimates
fitted_density <- density(tau_hat)
shift <- fitted_density$x[which.max(fitted_density$y)]

shifted_tau_hat <- tau_hat - shift

# hypothesis tests for each tau(whether they were significantly different from zero)
test_statistic <- as.vector(shifted_tau_hat / stderror_tau)
pval <- 1 - pnorm(abs(test_statistic))


prediction <- pval < 0.05
table(prediction, truth$effect.size !=1)


## if do not consider correlation between z scores

XtX <- t(design_matrix) %*% design_matrix
XtX_inv <- chol2inv(chol(XtX))
tau_hat_nocorr <- XtX_inv %*% t(design_matrix) %*% Z_score
resid_nocorr <- Z_score - design_matrix %*% tau_hat_nocorr
sigma2_hat_nocorr <- (t(resid_nocorr) %*% resid_nocorr) / nrow(resid_nocorr)
covar_tau_nocorr <- drop(sigma2_hat_nocorr) * XtX_inv
stderror_tau_nocorr <- sqrt(diag(covar_tau_nocorr))
fitted_density_nocorr <- density(tau_hat_nocorr)
shift_nocorr <- fitted_density_nocorr$x[which.max(fitted_density_nocorr$y)]

shifted_tau_hat_nocorr <- tau_hat_nocorr - shift_nocorr

test_statistic_nocorr <- as.vector(shifted_tau_hat_nocorr / stderror_tau_nocorr)
pval_nocorr <- 1 - pnorm(abs(test_statistic_nocorr))
prediction_nocorr <- pval_nocorr < 0.05
table(prediction_nocorr, truth$effect.size !=1)

## in constrast with DeSeq2

library(DESeq2)
count_data <- as.matrix(data_withDA$obs.abn) |> unname()

metadata <- data.frame(id = sprintf("sub%d", seq(1, 100)),
                       group = c(rep(0, 50), rep(1, 50)))
metadata$group <- as.factor(metadata$group)

dds <- DESeqDataSetFromMatrix(countData=count_data,
                              colData=metadata,
                              design=~group)

dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_prediction <- deseq_result$padj < 0.05
table(deseq_prediction, truth$effect.size !=1)
