library(ExperimentHub)
library(HMP16SData)
library(GUniFrac)
library(phyloseq)
library(dplyr)


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

