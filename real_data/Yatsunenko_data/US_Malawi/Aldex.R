
library(ALDEx2)
rm(list=ls())


otu_infant <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_OTU.csv",
                       row.names = 1)
# prevalence_infant <- rowMeans(otu_infant != 0)
metadata_infant <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_metadata.csv",
                            row.names=1)
metadata_infant$Malawi <- metadata_infant$country == "MA"
otu_infant_filtered <- otu_infant[-c(4, 10), ]

result <- aldex(as.matrix(otu_infant_filtered),
                metadata_infant$country, denom="all")

aldex_summary <- data.frame(Taxa_Name = rownames(otu_infant_filtered),
                            pval = result$wi.eBH,
                            DA = result$wi.eBH < 0.05)

write.csv(aldex_summary, "real_data/Yatsunenko_data/US_Malawi/US_MA_infant_aldex_result.csv",
          row.names = FALSE)


