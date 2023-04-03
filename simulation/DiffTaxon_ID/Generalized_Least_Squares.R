## ----setup, include=FALSE----------------------------------------------------------
rm(list=ls())
library(Matrix)
library(ggplot2)
library(DESeq2)
library(cowplot)


## ----------------------------------------------------------------------------------
folder <- "/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation"
data <- readRDS(file.path(folder, 'data', 'PoiGamma',
                          'simple_simulation_withDA2.rds'))
truth <- data$mean.eco.abn
observed_counts <- data$obs.abn
num_taxa <- nrow(truth) # 100 taxa
glm_result <- read.csv(file.path(folder, 'glmdisp_result',
                                 'glmdisp_simple_result_withDA2.csv'))


## ----------------------------------------------------------------------------------
Z_score <- glm_result$effect / glm_result$SE
cat(sprintf("The number of response variables are %d", length(Z_score)))


## ----------------------------------------------------------------------------------
# one row for each taxa pair
allpairs <- combn(seq(1, num_taxa), m=2)
num_pair <- ncol(allpairs)
design_matrix_rowid <- rep(seq(1, num_pair), 2)
design_matrix_colid <- c(allpairs[1,], allpairs[2,])
entry_vals <- rep(c(1, -1), each=num_pair)

design_matrix <- sparseMatrix(i=design_matrix_rowid, j=design_matrix_colid,
                              x=entry_vals)

dm_rank <- rankMatrix(design_matrix, method='qr')[[1]]
cat(sprintf("Design matrix has a rank of %d. \n", dm_rank))



## ----------------------------------------------------------------------------------
nonzero_proportion <- rowMeans(observed_counts != 0)
taxa_reffect_var <- 1/nonzero_proportion # random effect variance
residual_var <- rep(1/4, num_pair) # residual variance


# Woodbury inverse
woodbury_inverse <- function(diagA, U_mat, diagC, V_mat){
  # calculate (A + UCV)^-1, A and C being diagonal matrix
  A_inv <- sparseMatrix(i=1:length(diagA), j=1:length(diagA),
                        x=1/diagA) # feed in a vector, get a sparse matrix
  C_inv <- sparseMatrix(i=1:length(diagC), j=1:length(diagC),
                        x=1/diagC)

  C_VAU <- C_inv + V_mat %*% A_inv %*% U_mat
  C_VAU_inv <- chol2inv(chol(C_VAU))
  C_VAU_inv <- drop0(C_VAU_inv, tol=1e-12)
  output <- A_inv - A_inv %*% U_mat %*% C_VAU_inv %*% V_mat %*% A_inv
  output <- drop0(output, tol=1e-12)
  return(output)
}

# covariance matrix inverse(without the unknown sigma^2 factor)
V_inv <- woodbury_inverse(residual_var, design_matrix, taxa_reffect_var,
                          t(design_matrix))


## ----------------------------------------------------------------------------------
XtVX <- t(design_matrix) %*% V_inv %*% design_matrix
XtVX <- as.matrix(XtVX)

symmetric_pseudo_inverse <- function(symmat, rank){
  # Moore-Penrose inverse for positive semidefinite matrix with no full rank
  eigen_result <- eigen(symmat, symmetric = TRUE)
  eigen_values <- eigen_result$values
  eigen_vectors <- eigen_result$vectors
  inv_eig_values <- rep(0, length(eigen_values))
  inv_eig_values[1:rank] <- 1 / eigen_values[1:rank]
  inv_mat <- eigen_vectors %*% diag(inv_eig_values) %*% t(eigen_vectors)
  return(inv_mat)
}

XtVX_inv <- symmetric_pseudo_inverse(XtVX, rank=dm_rank)

## ----------------------------------------------------------------------------------
# estimate tau
tau_hat <- XtVX_inv %*% t(design_matrix) %*% V_inv %*% Z_score
tau_hat <- as.vector(tau_hat)
# residual
resid <- Z_score - design_matrix %*% tau_hat
# sigma^2 estimate (MLE)
sigma2_hat <- t(resid) %*% V_inv %*% resid / nrow(resid)
# estimated covariance of tau
covar_tau <- drop(sigma2_hat) * XtVX_inv
stderror_tau <- sqrt(diag(covar_tau))


## ---- fig.width=7, fig.height=3----------------------------------------------------
strength <- data.frame(tau_hat = tau_hat)
strength$logfoldchange <- -log(truth$effect.size)

histogram_plot <- ggplot(strength, aes(x=tau_hat)) +
  geom_histogram(color="black", fill="white", binwidth=1)+
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  xlab("Estimated Tau") + ylab("Frequency") + theme_bw()

scatter_plot <- ggplot(strength, aes(x=tau_hat, y=logfoldchange)) +
  geom_point(size=0.8, alpha=0.8) +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  xlab("Estimated Tau") +
  ylab("Log Fold Change of Absolute Abundance") +
  theme_bw()

plot_grid(histogram_plot, scatter_plot)


## ----------------------------------------------------------------------------------
shift <- median(tau_hat) # shift the estimated taus so that the median is zero
shifted_tau_hat <- tau_hat - shift


## ----------------------------------------------------------------------------------
test_statistic <- as.vector(shifted_tau_hat / stderror_tau)
pval <- pnorm(-abs(test_statistic))
adjusted_pval <- p.adjust(pval, method="BH")
performance <- truth
performance$prediction <- adjusted_pval < 0.05
contingency_table <- table(performance$prediction, truth$effect.size !=1)


## ----------------------------------------------------------------------------------
# OLS estimates
XtX <- t(design_matrix) %*% design_matrix
XtX_inv <- symmetric_pseudo_inverse(XtX, rank=dm_rank)
tau_hat_nocorr <- XtX_inv %*% t(design_matrix) %*% Z_score
resid_nocorr <- Z_score - design_matrix %*% tau_hat_nocorr
sigma2_hat_nocorr <- (t(resid_nocorr) %*% resid_nocorr) / nrow(resid_nocorr)
covar_tau_nocorr <- drop(sigma2_hat_nocorr) * XtX_inv
stderror_tau_nocorr <- sqrt(diag(covar_tau_nocorr))

# shift
shift_nocorr <- median(tau_hat_nocorr)
shifted_tau_hat_nocorr <- tau_hat_nocorr - shift_nocorr

# test
test_statistic_nocorr <- as.vector(shifted_tau_hat_nocorr / stderror_tau_nocorr)
pval_nocorr <- pnorm(-abs(test_statistic_nocorr))
pval_nocorr_adjusted <- p.adjust(pval_nocorr)
performance$prediction_nocorr <- pval_nocorr_adjusted < 0.05
contingency_table_nocorr <- table(performance$prediction_nocorr,
                                  truth$effect.size !=1)



## ---- message=FALSE----------------------------------------------------------------

count_data <- as.matrix(data$obs.abn) |> unname()
metadata <- data.frame(id = sprintf("sub%d", seq(1, 100)),
                       group = c(rep(0, 50), rep(1, 50)))
metadata$group <- as.factor(metadata$group)
dds <- DESeqDataSetFromMatrix(countData=count_data,
                              colData=metadata,
                              design=~group)
dds <- DESeq(dds)
deseq_result <- results(dds)
performance$deseq_prediction <- deseq_result$padj < 0.05
contingency_table_deseq <- table(performance$deseq_prediction,
                                 truth$effect.size !=1)



