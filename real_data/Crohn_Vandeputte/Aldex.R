
library(ALDEx2)
rm(list=ls())


otu_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_otu.csv",
                   row.names = 1)
otu_mat <- as.matrix(otu_df)
# prevalence <- rowMeans(otu_df != 0)
metadata_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_metadata.csv",
                        row.names=1)
metadata <- metadata_df$Crohn

result <- aldex(otu_mat, metadata, denom="iqlr")

aldex_summary <- data.frame(Taxa_Name = rownames(otu_mat),
                            pval = result$wi.eBH,
                            DA = result$wi.eBH < 0.05)

write.csv(aldex_summary, "real_data/Crohn_Vandeputte/Crohn_Aldex2_result.csv",
          row.names = FALSE)


