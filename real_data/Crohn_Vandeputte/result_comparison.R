library(ggplot2)
library(reshape2)
library(cowplot)
rm(list=ls())

otu_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_otu.csv",
                   row.names = 1)

metadata_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_metadata.csv",
                        row.names=1)

t_otu <- t(otu_df) |> as.data.frame()
t_otu <- data.frame(t_otu > 0)

t_otu$Crohn <- metadata_df$Crohn
t_otu$Individual <- rownames(t_otu)
t_otu$Crohn <- as.factor(t_otu$Crohn)
long_otu <- melt(t_otu, id.vars=c("Crohn", "Individual"),
                 variable.name = "Genus")
ggplot(long_otu, aes(Individual, Genus, fill=Crohn, alpha=value))+
  geom_tile(colour = "gray50") +theme_bw() + xlab("Individual") + ylab("Genus")+ guides(alpha="none")+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())

ancomplus_result <- read.csv("real_data/Crohn_Vandeputte/Crohn_DiffRatio_result.csv")

ancombc_result <- read.csv("real_data/Crohn_Vandeputte/Crohn_ANCOMBC_result.csv")
aldex_result <- read.csv("real_data/Crohn_Vandeputte/Crohn_Aldex2_result.csv")
dacomp_result <- read.csv("real_data/Crohn_Vandeputte/Crohn_DACOMP_result.csv")
dacomp_result$DA[is.na(dacomp_result$DA)] <- FALSE


ancomplus_result$tau_lower <- ancomplus_result$tau_hat - 1.96*ancomplus_result$stderror_tau
ancomplus_result$tau_higher <- ancomplus_result$tau_hat + 1.96*ancomplus_result$stderror_tau
forest_plot <- ggplot(data=ancomplus_result, aes(x=Taxa, y=tau_hat, ymin=tau_lower, ymax=tau_higher, color=DA)) +
  geom_pointrange() +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Genus") + ylab("Estimated Tau (95% CI)") +
  theme_bw() + theme(legend.position = "none",axis.text.y=element_blank())# use a white background

histogram_plot <- ggplot(data=ancomplus_result, aes(x=tau_hat))+
  geom_histogram(color="black", fill="white", binwidth=0.5) + xlab("Estimated Tau")+
  geom_vline(xintercept=0, lty=2)+theme_bw()

plot_grid(forest_plot, histogram_plot, nrow=1)


combined_results <- cbind(ancomplus_result$Taxa,
                          ancomplus_result$DA,
                          aldex_result$DA,
                          ancombc_result$DA,
                          dacomp_result$DA) |> as.data.frame()
colnames(combined_results) <- c("Genus", "DiffRatio", "Aldex2", "ANCOMBC",
                                "DACOMP")


long_combined_results <- melt(combined_results, id.vars = "Genus",
                              variable.name = "Method")

ggplot(long_combined_results, aes(Genus, Method, alpha=value))+
  geom_tile(colour = "gray50") +theme_bw() + xlab("Genus") + ylab("Method")+ guides(alpha="none")+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank())



