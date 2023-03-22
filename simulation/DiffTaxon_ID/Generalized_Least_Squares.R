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
num_taxa <- nrow(truth) # 100 taxa
glm_result <- read.csv(file.path(folder, 'glmdisp_result', 
                                 'glmdisp_simple_result_withDA2.csv'))


## ----------------------------------------------------------------------------------
Z_score <- glm_result$effect / glm_result$SE
Z_score <- c(0, Z_score)
cat(sprintf("The number of response variables are %d", length(Z_score)))


## ----------------------------------------------------------------------------------
# one row for each taxa pair
allpairs <- combn(seq(1, num_taxa), m=2)
design_matrix <- matrix(0, nrow=nrow(glm_result), ncol=num_taxa)
for (j in 1:ncol(allpairs)){
  design_matrix[j, allpairs[1, j]] <- 1
  design_matrix[j, allpairs[2, j]] <- -1
}
initial_rank <- rankMatrix(design_matrix)[[1]]
cat(sprintf("Initial design matrix has a rank of %d. \n", initial_rank))

# one extra row for theta_1
design_matrix <- rbind(c(1, rep(0, num_taxa - 1)), design_matrix)
final_rank <- rankMatrix(design_matrix)[[1]]
cat(sprintf("The final design matrix has a rank of %d.", final_rank))


## ----------------------------------------------------------------------------------
V_mat <- design_matrix %*% t(design_matrix)
V_mat <- V_mat + diag(1/4, nrow=nrow(V_mat), ncol=ncol(V_mat))
V_mat[1,1] <- 1 # the top left entry doesn't have epsilon^2


## ----------------------------------------------------------------------------------
if (file.exists(file.path(folder, 'DiffTaxon_ID', 'example_V_inv.rds'))){
  V_inv <- readRDS(file.path(folder, 'DiffTaxon_ID', 'example_V_inv.rds'))
} else{
  V_inv <- chol2inv(chol(V_mat))
}


## ----------------------------------------------------------------------------------
XtVX <- t(design_matrix) %*% V_inv %*% design_matrix
XtVX_inv <- chol2inv(chol(XtVX))


## ----------------------------------------------------------------------------------
# estimate tau
tau_hat <- XtVX_inv %*% t(design_matrix) %*% V_inv %*% Z_score
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
fitted_density <- density(tau_hat)
shift <- fitted_density$x[which.max(fitted_density$y)]
shifted_tau_hat <- tau_hat - shift


## ----------------------------------------------------------------------------------
test_statistic <- as.vector(shifted_tau_hat / stderror_tau)
pval <- pnorm(-abs(test_statistic))
performance <- truth
performance$prediction <- pval < 0.05
contingency_table <- table(performance$prediction, truth$effect.size !=1)


## ----------------------------------------------------------------------------------
# OLS estimates
XtX <- t(design_matrix) %*% design_matrix
XtX_inv <- chol2inv(chol(XtX))
tau_hat_nocorr <- XtX_inv %*% t(design_matrix) %*% Z_score
resid_nocorr <- Z_score - design_matrix %*% tau_hat_nocorr
sigma2_hat_nocorr <- (t(resid_nocorr) %*% resid_nocorr) / nrow(resid_nocorr)
covar_tau_nocorr <- drop(sigma2_hat_nocorr) * XtX_inv
stderror_tau_nocorr <- sqrt(diag(covar_tau_nocorr))

# shift
fitted_density_nocorr <- density(tau_hat_nocorr)
shift_nocorr <- fitted_density_nocorr$x[which.max(fitted_density_nocorr$y)]
shifted_tau_hat_nocorr <- tau_hat_nocorr - shift_nocorr

# test
test_statistic_nocorr <- as.vector(shifted_tau_hat_nocorr / stderror_tau_nocorr)
pval_nocorr <- pnorm(-abs(test_statistic_nocorr))
performance$prediction_nocorr <- pval_nocorr < 0.05
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



