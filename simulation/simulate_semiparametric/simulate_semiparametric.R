library(ExperimentHub)
library(HMP16SData)
library(GUniFrac)
library(dirmult)
library(phyloseq)
library(dplyr)
library(MASS)

remove(list=ls())

summary_stats <- V35() %>% table_one()
colnames(summary_stats) <- c("Visit_Number", "Sex", "Center",
                             "Body_Site", "Body_Subsite")
summary_stats_stool <- summary_stats %>% filter(Body_Subsite == "Stool")

# load stool data in variable region 3-5
V35_stool <-
  V35() %>%
  subset(select = HMP_BODY_SUBSITE == "Stool")

# limit to unique patients at visit one
V35_stool_phyloseq <-
  as_phyloseq(V35_stool)
V35_stool_visitone <- subset_samples(V35_stool_phyloseq, VISITNO==1)

stool_taxonomy_table <- tax_table(V35_stool_visitone) %>% as.data.frame()
stool_sample_table <- sample_data(V35_stool_visitone)
stool_count_table <- otu_table(V35_stool_visitone) %>% as.matrix()

# select the top 500 taxa in terms of prevalence
existence <- stool_count_table > 0
summary_existence <- rowMeans(existence) %>% sort(decreasing = TRUE)
taxa_included <- names(summary_existence[1:500])


template_mat <- stool_count_table[taxa_included, ]
template_data <- template_mat@.Data
row.names(template_data) <- sprintf('Taxon%d', seq(1, 500))
colnames(template_data) <- sprintf('Indv%d', seq(1, 191))


folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/data/semiparametric'


# simulate null data
## estimate dirichlet multinomial distribution
rdirichlet.m <- function (alpha) {
  Gam <- rgamma(length(alpha), shape = alpha)
  Gam / sum(Gam)
}


simMSeq_null <- function(nSam=400, nTaxa=300, depth.mu=10000, depth.theta=100, seed=1){
  # simulate count matrix with no DA taxa
  set.seed(seed+2020)
  nSeq <- rnegbin(nSam, mu = depth.mu , theta = depth.theta) # generating sequencing depths
  # generate count matrices from multinomial distributions
  otu.tab.sim <- sapply(1:nSam, function (i) rmultinom(1, nSeq[i], rdirichlet.m(alpha = rep(0.5, nTaxa))))
  covariates <- rep(0, nSam)
  covariates[sample(nSam, nSam/2)] <- 1
  # filter out taxa that have too many zeros in either group
  group0 <- otu.tab.sim[, covariates == 0]
  group1 <- otu.tab.sim[, covariates == 1]
  admissible_taxa <- (rowMeans(group0 == 0) < 0.8) & (rowMeans(group1 == 0) < 0.8)
  otu.tab.sim <- otu.tab.sim[admissible_taxa, ]
  row.names(otu.tab.sim) <- sprintf("Taxon%d", seq(1, sum(admissible_taxa)))
  return(list(otu.tab.sim=otu.tab.sim,  covariates=covariates))
}



for (j in 1:100){
  simulated_data <- simMSeq_null(seed=j)
  output_file <- sprintf('simulated_data_null_%d.rds', j)
  saveRDS(simulated_data, file.path(folder, output_file))
}



# simulate datasets with DAs
for (j in 1:100){
  print(j)
  set.seed(j)
  simulated_data <- SimulateMSeq(ref.otu.tab = template_data,
                                 nSam=100, nOTU=500, diff.otu.pct = 0.2,
                                 diff.otu.direct = "balanced",
                                 diff.otu.mode = "mix",
                                 covariate.type = "binary",
                                 grp.ratio = 1,
                                 covariate.eff.mean = 1,
                                 covariate.eff.sd = 0.4,
                                 confounder.type = "none",
                                 depth.mu = 10000,
                                 depth.theta = 10)
  output_file <- sprintf('simulated_data_%d.rds', j)
  saveRDS(simulated_data, file.path(folder, output_file))
}

