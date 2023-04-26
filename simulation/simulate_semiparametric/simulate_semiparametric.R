library(ExperimentHub)
library(HMP16SData)
library(GUniFrac)
library(dirmult)
library(phyloseq)
library(dplyr)
library(MASS)
library(gtools)

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


folder <- '/home/wangmk/UM/Research/MDAWG/DiffRatio/simulation/'
source(file.path(folder, 'simulate_semiparametric', 'utils.R'))


cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)



simulated_datasets <- foreach(j=1:100, .inorder=FALSE,
                              .packages=c("gtools", "dirmult", "MASS")) %dopar% {

  simulated_data <- SimulateCount(template_data, nSample=100, diff_prop=0.2,
                                  covariate_type = "binary", grp_ratio=1, zinf_prop=0.1,
                                  covariate_eff_mean = 1.5, covariate_eff_sd = 0.2,
                                  depth_mu = 5000, depth_theta = 2, seed=j)


  return(simulated_data)
}


for (j in 1:100){
  filename <- sprintf("simulated_data_%d.rds", j)
  saveRDS(simulated_datasets[[j]], file.path(folder, 'data', 'semiparametric', filename))
}

