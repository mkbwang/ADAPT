library(dplyr)
library(ggplot2)
library(lme4)
library(zoib)


rm(list=ls())
data_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/RMD/CAARS_data'
model_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/RMD/CAARS_Model_Summary'
plot_folder = '/home/wangmk/UM/Research/MDAWG/Differential_Ratio/DiffRatio/RMD/CAARS_plots'

# load the count data and the metadata
count_order_level <- read.csv(file.path(data_folder, 'filtered_count_order_level.csv'),
                             row.names =1)
metadata <- read.csv(file.path(data_folder, 'filtered_metadata_order_level.csv'))

# plot all the pairs of taxa
taxa_names <- colnames(count_order_level)
taxa_pairs <- combn(taxa_names, 2) %>% t()

biplot_viz <- function(tpair){
  subset_data <- cbind(metadata$asthma, count_order_level[, tpair]) %>%
    as.data.frame()
  colnames(subset_data) <- c("Asthma", tpair)
  subset_data$ratios = subset_data[, 2] / (subset_data[, 2] + subset_data[, 3])
  subset_data$Asthma <- as.factor(subset_data$Asthma)

  plot1_name <- sprintf('%s_%s_biplot.pdf', tpair[1], tpair[2])
  count_plot <- ggplot(subset_data, aes_string(x=tpair[1], y=tpair[2], color='Asthma')) +
    geom_jitter(alpha=0.6) +
    theme_bw()
  ggsave(file.path(plot_folder, 'biplots', plot1_name),
         plot=count_plot,
         width=12, height=8, units="cm")

  plot2_name <- sprintf('%s_%s_Ratios.pdf', tpair[1], tpair[2])
  ratio_plot <- ggplot(subset_data, aes(x=ratios, color=Asthma)) +
    geom_histogram(fill="white")+
    xlab(sprintf('%s/(%s+%s)', tpair[1], tpair[1], tpair[2])) + theme_bw()
  ggsave(file.path(plot_folder, 'ratios', plot2_name),
         plot=ratio_plot,
         width=12, height=8, units="cm")
}

for (i in 1:nrow(taxa_pairs)){
  biplot_viz(taxa_pairs[i, ])
}

# an example pair that exposes the shortcoming of ANCOM and GLMM
example_pair1 <- c("o__Campylobacterales", "o__Synergistales")

example_data1 <- cbind(metadata$asthma, count_order_level[, example_pair1]) %>%
  as.data.frame()
# plot
colnames(example_data1) <- c("Asthma", example_pair1)
example_data1$Asthma <- as.factor(example_data1$Asthma)
ggplot(example_data1, aes_string(x=example_pair1[1], y=example_pair1[2], color="Asthma"))+
  geom_jitter(alpha=0.6) +
  theme_bw()
example_data1$ratios = example_data1[, 2] / (example_data1[, 2] + example_data1[, 3])
ggplot(example_data1, aes(x=ratios, color=Asthma)) +
  geom_histogram(fill="white")+
  xlab(sprintf('%s/(%s+%s)', example_pair1[1], example_pair1[1], example_pair1[2])) + theme_bw()

# wilcoxon test
example_data1$logratio <- log(example_data1[, 2] + 1) - log(example_data1[, 3] + 1)
wilx_test1 <- wilcox.test(logratio~Asthma, data=example_data1)

# GLMM
example_data1_copy <- example_data1
example_data1_copy$Asthma <- as.numeric(example_data1_copy$Asthma)-1
example_data1_copy$ID <- seq(1, nrow(example_data1_copy))
glmm1 <- glmer(cbind(o__Campylobacterales, o__Synergistales) ~ Asthma + (1|ID),
               data=example_data1_copy, family="binomial", nAGQ = 50,
               control=lme4::glmerControl(optimizer='bobyqa',
                                          optCtrl=list(maxfun=2e5)))

# zoib
if (file.exists(file.path(model_folder, 'oib_model.rds'))){
  oibmodel <- readRDS(file.path(model_folder, 'oib_model.rds'))
}else{
  oibmodel <- zoib(ratios ~ Asthma|Asthma|Asthma,
                   data=example_data1_copy, random=0, zero.inflation = FALSE,
                   one.inflation = TRUE, joint=FALSE,
                   n.iter=7000, n.thin=15, n.burn=1000)
  saveRDS(oibmodel, file=file.path(model_folder, 'oib_model.rds'))
}



coefficients <- oibmodel$coeff
summary(coefficients)


# example_pair2 <- c("o__Campylobacterales", "o__Flavobacteriales")
# example_pair3 <- c("o__Cardiobacteriales", "o__Saccharimonadales")


example_data2 <- cbind(metadata$asthma, count_order_level[, example_pair2]) %>%
  as.data.frame()
colnames(example_data2) <- c("Asthma", example_pair2)
example_data2$Asthma <- as.factor(example_data2$Asthma)
ggplot(example_data2, aes_string(x=example_pair2[1], y=example_pair2[2], color="Asthma"))+
  geom_jitter(alpha=0.6) +
  theme_bw()
example_data2$ratios = example_data2[, 2] / (example_data2[, 2] + example_data2[, 3])
ggplot(example_data2, aes(x=ratios, color=Asthma)) +
  geom_histogram(fill="white")+
  xlab(sprintf('%s/(%s+%s)', example_pair2[1], example_pair2[1], example_pair2[2])) + theme_bw()


example_data3 <- cbind(metadata$asthma, count_order_level[, example_pair3]) %>%
  as.data.frame()
colnames(example_data3) <- c("Asthma", example_pair3)
example_data3$Asthma <- as.factor(example_data3$Asthma)
example_data3$ratios = example_data3[, 2] / (example_data3[, 2] + example_data3[, 3])
ggplot(example_data3, aes_string(x=example_pair3[1], y=example_pair3[2], color="Asthma"))+
  geom_jitter(alpha=0.6) +
  theme_bw()
ggplot(example_data3, aes(x=ratios, color=Asthma)) +
  geom_histogram(fill="white")+
  xlab(sprintf('%s/(%s+%s)', example_pair3[1], example_pair3[1], example_pair3[2])) + theme_bw()

