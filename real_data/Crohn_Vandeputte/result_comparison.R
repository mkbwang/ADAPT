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

ancomplus_result_df <- read.csv("real_data/Crohn_Vandeputte/Crohn_DiffRatio_result.csv")
ancomplus_decision <- readRDS("real_data/Crohn_Vandeputte/Crohn_DiffRatio_result.rds")

ancombc_result <- readRDS("real_data/Crohn_Vandeputte/Crohn_ANCOMBC_result.rds")
aldex_result <- read.csv("real_data/Crohn_Vandeputte/Crohn_Aldex2_result.csv")
dacomp_result <- read.csv("real_data/Crohn_Vandeputte/Crohn_DACOMP_result.csv")
dacomp_result$DA[is.na(dacomp_result$DA)] <- FALSE


ancomplus_result_df$tau_lower <- ancomplus_result_df$tau_hat - 1.96*ancomplus_result_df$sderror
ancomplus_result_df$tau_higher <- ancomplus_result_df$tau_hat + 1.96*ancomplus_result_df$sderror
forest_plot <- ggplot(data=ancomplus_result_df, aes(x=Taxa_name, y=tau_hat, ymin=tau_lower, ymax=tau_higher)) +
  geom_pointrange(size=0.2) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Genus") + ylab("Estimated Tau (95% CI)") +
  theme_bw() + theme(legend.position = "none",axis.text.y=element_blank())# use a white background

histogram_plot <- ggplot(data=ancomplus_result_df, aes(x=tau_hat))+
  geom_histogram(color="black", fill="white", binwidth=0.5) + xlab("Estimated Tau")+
  geom_vline(xintercept=0, lty=2)+theme_bw()

plot_grid(forest_plot, histogram_plot, nrow=1)


combined_results <- cbind(ancomplus_result_df$Taxa_name,
                          ancomplus_result_df$Taxa_name %in% ancomplus_decision$difftaxa,
                          aldex_result$DA,
                          ancomplus_result_df$Taxa_name %in% ancombc_result,
                          dacomp_result$DA) |> as.data.frame()
colnames(combined_results) <- c("Genus", "DiffRatio", "Aldex2", "ANCOMBC",
                                "DACOMP")
combined_results$DiffRatio <- combined_results$DiffRatio == "TRUE"
combined_results$Aldex2 <- combined_results$Aldex2 == "TRUE"
combined_results$ANCOMBC <- combined_results$ANCOMBC == "TRUE"
combined_results$DACOMP <- combined_results$DACOMP == "TRUE"

DiffRatio_set <- combined_results$Genus[combined_results$DiffRatio]
Aldex2_set <- combined_results$Genus[combined_results$Aldex2]
ANCOMBC_set <- combined_results$Genus[combined_results$ANCOMBC]
DACOMP_set <- combined_results$Genus[combined_results$DACOMP]

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(DiffRatio_set, Aldex2_set, ANCOMBC_set, DACOMP_set),
  category.names = c("DiffRatio" , "Aldex2" , "ANCOMBC", "DACOMP"),
  filename = 'DAtaxa_overlap.pdf',
  output=TRUE
)

long_combined_results <- melt(combined_results, id.vars = "Genus",
                              variable.name = "Method")

ggplot(long_combined_results, aes(Genus, Method, alpha=value))+
  geom_tile(colour = "gray50") +theme_bw() + xlab("Genus") + ylab("Method")+ guides(alpha="none")+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank())



