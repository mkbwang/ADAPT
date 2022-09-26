rm(list=ls())
folder = '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation'
library(dplyr)


data <- readRDS(file.path(folder, 'simulated_data.rds'))
abn_info <- data$mean.eco.abn

taxa_names <- row.names(abn_info)
taxa_pairs <- combn(taxa_names, 2) %>% t() %>% as.data.frame()
colnames(taxa_pairs) <- c("T1", "T2")

# collect the truth
truth <- taxa_pairs
truth$diffratio <- TRUE
for (j in 1:nrow(truth)){
  t1 <- truth$T1[j]
  t2 <- truth$T2[j]
  truth$diffratio[j] <- abn_info[t1, 'effect.size'] != abn_info[t2, 'effect.size']
}

# load model results
ANCOM_result <- readRDS(file.path(folder, 'ancom_result.rds'))
GLMM_result <- readRDS(file.path(folder, 'glmm_result.rds'))
ZOIB_result <- readRDS(file.path(folder, 'ZOIB_result.rds'))

combination <- cbind(truth, ANCOM_result$pval, GLMM_result$pval, ZOIB_result$zinfpval,
                     ZOIB_result$oinfpval, ZOIB_result$betapval)

colnames(combination) <- c("T1", "T2", "diffratio", "ANCOM_pval",
                           "GLMM_pval", 'ZOIB_zinf_pval', 'ZOIB_oinf_pval',
                           'ZOIB_beta_pval')


combination_replicate <- combination

combination_replicate$ANCOM_pval <- p.adjust(combination_replicate$ANCOM_pval,
                                             method='BH')
combination_replicate$GLMM_pval <- p.adjust(combination_replicate$GLMM_pval,
                                            method='BH')
combination_replicate$ZOIB_zinf_pval <- p.adjust(combination_replicate$ZOIB_zinf_pval,
                                            method='BH')
combination_replicate$ZOIB_oinf_pval <- p.adjust(combination_replicate$ZOIB_oinf_pval,
                                                 method='BH')
combination_replicate$ZOIB_beta_pval <- p.adjust(combination_replicate$ZOIB_beta_pval,
                                                 method='BH')


truth <- combination_replicate$diffratio
ancom_decision <- combination_replicate$ANCOM_pval <0.05
GLMM_decision <- combination_replicate$GLMM_pval <0.05

zinf_decision <- combination_replicate$ZOIB_zinf_pval < 0.05
zinf_decision[is.na(zinf_decision)] <- FALSE
oinf_decision <- combination_replicate$ZOIB_oinf_pval < 0.05
oinf_decision[is.na(oinf_decision)] <- FALSE
beta_decision <- combination_replicate$ZOIB_beta_pval < 0.05

ZOIB_decision <- zinf_decision | oinf_decision | beta_decision

table(truth, ancom_decision)
table(truth, GLMM_decision)
table(truth, ZOIB_decision)
table(GLMM_decision, ZOIB_decision)
table(ancom_decision, GLMM_decision)
table(ancom_decision, ZOIB_decision)



combination$ANCOM_decision <- ancom_decision
combination$GLMM_decision <- GLMM_decision
combination$ZOIB_decision <- ZOIB_decision


pairs1 <- combination %>% filter(!diffratio & !GLMM_decision & ZOIB_decision)
pairs2 <- combination %>% filter(diffratio & !GLMM_decision & ZOIB_decision)

sample_mat <- data$obs.abn


## example 1
example_1 <- sample_mat[c("taxon2", "taxon36"), ] %>% t() %>% as.data.frame()
example_1$grp <- data$grp - 1
example_1$prop <- example_1$taxon2/(example_1$taxon2 + example_1$taxon36)
example_1_copy <- example_1
example_1_copy$grp <- as.factor(example_1_copy$grp)


library(ggplot2)
library(cowplot)
library(betareg)
library(lme4)

scatterplot1 <- ggplot(example_1_copy, aes(x=taxon36, y=taxon2, color=grp)) +
  geom_jitter(size=2, alpha=0.4, width=0.2, height=0.2) +
  theme_bw() + theme(text = element_text(size = 14))
histogram1 <- ggplot(example_1_copy , aes(x=prop, fill = grp))+
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity') +
  theme_bw() + xlab("Proportion of taxon 2") + xlim(-0.05, 1.05)+
  theme(text = element_text(size = 14))

combined_plot1 <- plot_grid(scatterplot1, histogram1, align='h')

example_1$tot <- example_1$taxon2 + example_1$taxon36

example_1$ID <- row.names(example_1)
example_1_filtered <- example_1 %>% filter(taxon2 > 0 & taxon36 > 0)
example_1_filtered$weights <- example_1_filtered$tot/sum(example_1_filtered$tot) * nrow(example_1_filtered)

beta_example1 <- betareg(prop ~ grp|grp, weights=weights, data=example_1_filtered)
glmm_example1 <- glmer(cbind(taxon2, taxon36) ~ grp + (1|ID), family="binomial", data=example_1)




## second example

example_2 <- sample_mat[c("taxon152", "taxon161"), ] %>% t() %>% as.data.frame()
example_2$grp <- data$grp - 1
example_2$prop <- example_2$taxon152/(example_2$taxon152 + example_2$taxon161)
example_2_copy <- example_2
example_2_copy$grp <- as.factor(example_2_copy$grp)

scatterplot2 <- ggplot(example_2_copy, aes(x=taxon161, y=taxon152, color=grp)) +
  geom_jitter(size=2, alpha=0.4, width=0.2, height=0.2) +
  theme_bw() + theme(text = element_text(size = 14))
histogram2 <- ggplot(example_2_copy , aes(x=prop, fill = grp))+
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity') +
  theme_bw() + xlab("Proportion of taxon 152") + xlim(-0.05, 1.05)+
  theme(text = element_text(size = 14))

combined_plot2 <- plot_grid(scatterplot2, histogram2, align='h')

example_2$tot <- example_2$taxon152 + example_2$taxon161
example_2$ID <- row.names(example_2)
example_2_filtered <- example_2 %>% filter(taxon152 > 0 & taxon161 > 0)
example_2_filtered$weights <- example_2_filtered$tot/sum(example_2_filtered$tot) * nrow(example_2_filtered)

beta_example2 <- betareg(prop ~ grp|grp, weights=weights, data=example_2_filtered)
glmm_example2 <- glmer(cbind(taxon152, taxon161) ~ grp + (1|ID), family="binomial", data=example_2,
                       nAGQ=20)


