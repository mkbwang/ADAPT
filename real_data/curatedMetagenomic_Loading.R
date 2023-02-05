# curatedMetagenomicData trial

library(ExperimentHub)
library(HMP16SData)
library(curatedMetagenomicData)
library(phyloseq)
library(dplyr)

summary_stats <- V35() %>% table_one()

location_counts <- summary_stats %>% group_by(`HMP Body Subsite`) %>% summarise(count=n())

V35_stool <-
  V35() %>%
  subset(select = HMP_BODY_SUBSITE == "Stool")


V35_stool_phyloseq <-
  as_phyloseq(V35_stool)

stool_taxonomy_table <- tax_table(V35_stool_phyloseq) %>% as.data.frame()

stool_count_table <- otu_table(V35_stool_phyloseq) %>% as.matrix()
existence <- stool_count_table > 0
summary_existence <- rowMeans(existence)


diet_microbiome <- sampleMetadata %>% filter(study_name == 'AsnicarF_2021')
diet_microbiome_abundance <- returnSamples(diet_microbiome, dataType = 'relative_abundance',
                                           counts=TRUE)
diet_microbiome_abundance <- mia::makePhyloseqFromTreeSummarizedExperiment(diet_microbiome_abundance,
                                                                           abund_values = 'relative_abundance')

diet_microbiome_count <- otu_table(diet_microbiome_abundance)

diet_microbiome_metadata <- sample_data(diet_microbiome_abundance)


Cameroon_metadata <- sampleMetadata %>% filter(study_name == 'RubelMA_2020')

Cameroon_abundance <- returnSamples(Cameroon_metadata, dataType = 'relative_abundance',
                                    counts=TRUE)
Cameroon_abundance <- mia::makePhyloseqFromTreeSummarizedExperiment(Cameroon_abundance,
                                                                    abund_values = "relative_abundance")

Cameroon_taxonomy <- tax_table(Cameroon_abundance) %>% as.data.frame()


Irish_metadata <- sampleMetadata %>% filter(study_name == 'KeohaneDM_2020')

Irish_abundance <- returnSamples(Irish_metadata, dataType = 'relative_abundance',
                                    counts=TRUE)
Irish_abundance <- mia::makePhyloseqFromTreeSummarizedExperiment(Irish_abundance,
                                                                    abund_values = "relative_abundance")

Irish_taxonomy <- tax_table(Irish_abundance) %>% as.data.frame()


HMP_2012_metadata <- sampleMetadata %>% filter(study_name == 'HMP_2012')
HMP_2012_abundance <- returnSamples(HMP_2012_metadata, dataType = 'relative_abundance',
                                    counts=TRUE)
HMP_2012_abundance <- mia::makePhyloseqFromTreeSummarizedExperiment(HMP_2012_abundance,
                                                                    abund_values = "relative_abundance")


HMP_taxonomy <- tax_table(HMP_2012_abundance) %>% as.data.frame()


IBD <- sampleMetadata %>% filter(study_name == 'VilaAV_2018')

lifelinesdeep_metadata <- sampleMetadata %>% filter(study_name == 'LifeLinesDeep_2016')
lifelinesdeep_abundance <- returnSamples(lifelinesdeep_metadata, dataType = 'relative_abundance',
                                    counts=TRUE)
lifelinesdeep_abundance <- mia::makePhyloseqFromTreeSummarizedExperiment(lifelinesdeep_abundance,
                                                                    abund_values = "relative_abundance")

sample_metadata <- sam_data(lifelinesdeep_abundance)

obesity <- sampleMetadata %>% filter(study_name == 'LeChatelierE_2013')

obesity_abundance <- returnSamples(obesity, dataType = 'relative_abundance',
                                   counts=TRUE)
obesity_abundance <- mia::makePhyloseqFromTreeSummarizedExperiment(obesity_abundance,
                                                                         abund_values = "relative_abundance")
obesity_count_table <- otu_table(obesity_abundance)
obesity_taxa <- tax_table(obesity_abundance)


