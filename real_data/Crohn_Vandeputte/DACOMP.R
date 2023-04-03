library(ALDEx2)
rm(list=ls())


otu_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_otu.csv",
                   row.names = 1)
otu_mat <- as.matrix(otu_df)
# prevalence <- rowMeans(otu_df != 0)
metadata_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_metadata.csv",
                        row.names=1)
metadata <- metadata_df$Crohn


# find reference taxa
result.selected.references = dacomp.select_references(
  X = t(otu_mat), verbose = F)

# DA test for the rest of taxa
result.test = dacomp.test(X = t(otu_mat), #counts data
                          y = metadata, #phenotype in y argument
                          # obtained from dacomp.select_references(...):
                          ind_reference_taxa = result.selected.references,
                          test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                          verbose = F)

pval_adjusted <- result.test$p.values.test.adjusted
decision <- result.test$dsfdr_rejected

DACOMP_summary <- data.frame(Taxa_Name = rownames(otu_mat),
                             Pval_adjust = pval_adjusted,
                             DA = pval_adjusted < 0.05)


write.csv(DACOMP_summary, "real_data/Crohn_Vandeputte/Crohn_DACOMP_result.csv",
          row.names = FALSE)



