library(ggplot2)
library(reshape2)

rm(list=ls())
otu_infant <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_OTU.csv",
                       row.names = 1)
otu_infant_filtered <- otu_infant[-c(4, 10), ]
t_otu_infant <- t(otu_infant_filtered) |> as.data.frame()
t_otu_infant <- data.frame(t_otu_infant > 0)


metadata_infant <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_metadata.csv",
                            row.names=1)
t_otu_infant$Country <- metadata_infant$country
t_otu_infant$Infant <- rownames(t_otu_infant)
long_otu_infant <- melt(t_otu_infant, id.vars=c("Country", "Infant"),
                        variable.name = "Phylum")
long_otu_infant$value <- as.numeric(long_otu_infant$value)
ggplot(long_otu_infant, aes(Infant, Phylum, fill=Country, alpha=value))+
  geom_tile(colour = "gray50") +theme_bw() + xlab("Infant") + ylab("Phylum")+ guides(alpha="none")+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())


ancomplus_result <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_DiffRatio_result.csv")

ancombc_result <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_ANCOMBC_result.csv")
aldex_result <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_aldex_result.csv")
dacomp_result <- read.csv("real_data/Yatsunenko_data/US_Malawi/US_MA_infant_DACOMP_result.csv")
dacomp_result$DA[is.na(dacomp_result$DA)] <- FALSE


ancomplus_result$tau_lower <- ancomplus_result$tau_hat - 1.96*ancomplus_result$stderror_tau
ancomplus_result$tau_higher <- ancomplus_result$tau_hat + 1.96*ancomplus_result$stderror_tau
ggplot(data=ancomplus_result, aes(x=Taxa, y=tau_hat, ymin=tau_lower, ymax=tau_higher, color=DA)) +
  geom_pointrange() +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Phylum") + ylab("Estimated Tau (95% CI)") +
  theme_bw() +theme(legend.position = "None")  # use a white background



combined_results <- cbind(ancomplus_result$Taxa,
                          ancomplus_result$DA,
                          aldex_result$DA,
                          ancombc_result$DA,
                          dacomp_result$DA) |> as.data.frame()
colnames(combined_results) <- c("Phylum", "DiffRatio", "Aldex2", "ANCOMBC",
                                "DACOMP")


long_combined_results <- melt(combined_results, id.vars = "Phylum",
                              variable.name = "Method")

ggplot(long_combined_results, aes(Phylum, Method, alpha=value))+
  geom_tile(colour = "gray50") +theme_bw() + xlab("Phylum") + ylab("Method")+ guides(alpha="none")+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))



