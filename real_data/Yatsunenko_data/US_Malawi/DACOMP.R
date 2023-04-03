library(dacomp)
rm(list=ls())


otu_infant <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_OTU.csv",
                       row.names = 1)
# prevalence_infant <- rowMeans(otu_infant != 0)
metadata_infant <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_metadata.csv",
                            row.names=1)
metadata_infant$Malawi <- metadata_infant$country == "MA"
otu_infant_filtered <- otu_infant[-c(4, 10), ]
otu_infant_filtered_mat <- as.matrix(otu_infant_filtered)


# find reference taxa
result.selected.references = dacomp.select_references(
  X = t(otu_infant_filtered_mat), verbose = F)

# DA test for the rest of taxa
result.test = dacomp.test(X = t(otu_infant_filtered_mat), #counts data
                          y = metadata_infant$country, #phenotype in y argument
                          # obtained from dacomp.select_references(...):
                          ind_reference_taxa = result.selected.references,
                          test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                          verbose = F)

pval_adjusted <- result.test$p.values.test.adjusted
decision <- result.test$dsfdr_rejected

DACOMP_summary <- data.frame(Taxa_Name = rownames(otu_infant_filtered_mat),
                             Pval_adjust = pval_adjusted,
                             DA = pval_adjusted < 0.05)


write.csv(DACOMP_summary, "real_data/Yatsunenko_data/US_Malawi/US_MA_infant_DACOMP_result.csv",
          row.names = FALSE)



