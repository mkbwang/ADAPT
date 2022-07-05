
# load data
source('RMD/preprocess.R')
folder = 'RMD/CAARS_data'
load(file.path(file.path(folder, 'CAARS_processed_GENUS.Rdata')))
processed_data <- preprocess(CAARS.data.genus, 'SAMPLE_ID', "asthma")


filtered_count <- processed_data$feature_table
filtered_metadata <- processed_data$meta_data
struc_zero <- processed_data$structure_zeros

# check out the fisher test result
fisher_test_result <- read.csv('RMD/CAARS_Model_Summary/fisher_test.csv')
significant_pairs_ID <- which(fisher_test_result$adjustedp < 0.05)

# subset the significant pairs
significant_pairs <- fisher_test_result[significant_pairs_ID, ]

# pick one example to take a closer look
examplepair <- significant_pairs[6, c(1,2)] %>% unlist()
example_counts <- filtered_count[examplepair, ]

asthma_counts <- example_counts[, filtered_metadata$asthma == 1]
nonasthma_counts <- example_counts[, filtered_metadata$asthma == 0]

asthma_dbpos <- sum((asthma_counts[1,] > 0) & (asthma_counts[2,] > 0), na.rm=TRUE)
asthma_dbneg <- sum((asthma_counts[1,] == 0) & (asthma_counts[2,] == 0), na.rm=TRUE)
asthma_t1post2neg <- sum((asthma_counts[1,] > 0) & (asthma_counts[2,] == 0), na.rm=TRUE)
asthma_t1negt2pos <- sum((asthma_counts[1,] == 0) & (asthma_counts[2,] > 0), na.rm=TRUE)

control_dbpos <- sum((nonasthma_counts[1, ] > 0) & (nonasthma_counts[2, ] > 0), na.rm=TRUE)
control_dbneg <- sum((nonasthma_counts[1, ] == 0) & (nonasthma_counts[2, ] == 0), na.rm=TRUE)
control_t1post2neg <- sum((nonasthma_counts[1, ] > 0) & (nonasthma_counts[2, ] == 0), na.rm=TRUE)
control_t1negt2pos <- sum((nonasthma_counts[1, ] == 0) & (nonasthma_counts[2, ] > 0), na.rm=TRUE)

contingency_table <- matrix(c(asthma_dbpos, asthma_dbneg, asthma_t1post2neg, asthma_t1negt2pos,
                              control_dbpos, control_dbneg, control_t1post2neg, control_t1negt2pos),
                            nrow=4)


# GLMM

example_counts <- t(example_counts) %>% as.data.frame()
example_counts$asthma <- asthma
example_counts$ID <- row.names(example_counts)
example_counts <- example_counts %>% filter(g__Tannerella > 0 | g__Lentimicrobium > 0)

## directly use glmer function
library(lme4)
outcome1 <- glmer(cbind(g__Tannerella, g__Lentimicrobium) ~ asthma + (1|ID), data=example_counts, nAGQ=1,
                  family="binomial")


glm_result <- glm(cbind(g__Tannerella, g__Lentimicrobium) ~ asthma , data=example_counts, family="binomial")
X <- cbind(1, example_counts$asthma)
n <- example_counts$g__Tannerella + example_counts$g__Lentimicrobium
beta <- glm_result$coefficients %>% unname()
b <- rep(0, length(example_counts$asthma))
eta <- X %*% beta + b
pi <- 1/(1+exp(-eta)) %>% as.vector()
V <- diag(1/(n*pi*(1-pi)))
u <- eta + V %*% (example_counts$g__Tannerella - n * pi)

example_counts$link <- u
example_counts$glmvar <- 1/(n*pi*(1-pi))
example_counts$sqglmvar <- sqrt(1/(n*pi*(1-pi)))


## TODO: needs to figure out if nlme could work, I had better write out the function myself
library(nlme)

updated_model <- gls(link ~ asthma , data=example_counts,
                     weights=varConstProp(form=~sqglmvar))

