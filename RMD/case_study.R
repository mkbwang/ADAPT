library(lme4)

# case study

rm(list=ls())

data <- read.csv('RMD/CAARS_data/Lentimicrobium_Kingella.csv')
data_replicate <- data


## for ancom, we need to impute zeros first
lenti <- data$g__Lentimicrobium
lenti[lenti==0] <- 1
kingella <- data$g__Kingella
kingella[kingella==0] <- 1
data_replicate$g__Lentimicrobium <- lenti
data_replicate$g__Kingella <- kingella

data_replicate$logratio <- log(data_replicate$g__Lentimicrobium) - log(data_replicate$g__Kingella)

## Wilcoxon rank-sum test
wilx_test <- wilcox.test(logratio~asthma, data=data_replicate)
wilx_pval <- wilx_test$p.value

## T test
t_test <- t.test(logratio ~ asthma, data=data_replicate)
t_pval <- t_test$p.value

## Logistic Regression
glm_result <- glm(cbind(g__Lentimicrobium, g__Kingella) ~ asthma, data=data,
                  family = binomial(link = "logit"))
### hand calculation
asthma_patients <- data %>% filter(asthma == 1)
healthy_patients <- data %>% filter(asthma == 0)

healthy_lenti_count <- sum(healthy_patients$g__Lentimicrobium)
healthy_kingella_count <- sum(healthy_patients$g__Kingella)
asthma_lenti_count <- sum(asthma_patients$g__Lentimicrobium)
asthma_kingella_count <- sum(asthma_patients$g__Kingella)

glm_intercept <- log(healthy_lenti_count/
                       healthy_kingella_count)
sd_glm_intercept <- sqrt(1/healthy_lenti_count + 1/healthy_kingella_count)

glm_asthma_effect <- log(asthma_lenti_count/
                           asthma_kingella_count) -
  glm_intercept
sd_glm_asthma <- sqrt(1/healthy_lenti_count + 1/healthy_kingella_count+
                        1/asthma_lenti_count + 1/asthma_kingella_count)


## GLMM
glmm_result <- glmer(cbind(g__Lentimicrobium, g__Kingella) ~ asthma + (1|SAMPLE_ID), data=data,
                     family = binomial(link = "logit"), nAGQ=0L) # nAGQ is a big issue
