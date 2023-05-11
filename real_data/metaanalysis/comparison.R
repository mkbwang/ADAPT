
rm(list=ls())

baxter_performance <- read.csv('real_data/metaanalysis/baxter_performance.csv',
                               row.names=1)
baxter_polda_DA_rawp <- baxter_performance$Genus[baxter_performance$polda_rawp < 0.05] |> 
  na.omit()
baxter_polda_DA_adjustedp <- baxter_performance$Genus[baxter_performance$polda_adjustedp < 0.05] |> 
  na.omit()
baxter_polda_DA_rawp <- setdiff(baxter_polda_DA_rawp, baxter_polda_DA_adjustedp)
baxter_ancombc_DA_rawp <- baxter_performance$Genus[baxter_performance$ancombc_rawp < 0.05]
baxter_ancombc_DA_adjustedp <- baxter_performance$Genus[baxter_performance$ancombc_adjustedp < 0.05]
baxter_ancombc_DA_rawp <- setdiff(baxter_ancombc_DA_rawp, baxter_ancombc_DA_adjustedp)

baxter_aldex_DA_rawp <- baxter_performance$Genus[baxter_performance$aldex_rawp < 0.05]
baxter_aldex_DA_adjustedp <- baxter_performance$Genus[baxter_performance$aldex_adjustedp < 0.05]
baxter_aldex_DA_rawp <- setdiff(baxter_aldex_DA_rawp, baxter_aldex_DA_adjustedp)

baxter_result <- data.frame(POLDA_adjusted = rep('', 40),
                            POLDA_raw = rep('', 40),
                            ANCOMBC_adjusted = rep('', 40),
                            ANCOMBC_raw = rep('', 40),
                            ALDEX_adjusted = rep('', 40),
                            ALDEX_raw = rep('', 40))
baxter_result$POLDA_adjusted[1:length(baxter_polda_DA_adjustedp)] <- baxter_polda_DA_adjustedp
baxter_result$POLDA_raw[1:length(baxter_polda_DA_rawp)] <- baxter_polda_DA_rawp
baxter_result$ANCOMBC_adjusted[1:length(baxter_ancombc_DA_adjustedp)] <- baxter_ancombc_DA_adjustedp
baxter_result$ANCOMBC_raw[1:length(baxter_ancombc_DA_rawp)] <- baxter_ancombc_DA_rawp
baxter_result$ALDEX_raw[1:length(baxter_aldex_DA_rawp)] <- baxter_aldex_DA_rawp
baxter_result$ALDEX_adjusted[1:length(baxter_aldex_DA_adjustedp)] <- baxter_aldex_DA_adjustedp




zeller_performance <- read.csv('real_data/metaanalysis/zeller_performance.csv',
                               row.names=1)
zeller_polda_DA_rawp <- zeller_performance$Genus[zeller_performance$polda_rawp < 0.05] |> 
  na.omit()
zeller_polda_DA_adjustedp <- zeller_performance$Genus[zeller_performance$polda_adjustedp < 0.05] |> 
  na.omit()
zeller_polda_DA_rawp <- setdiff(zeller_polda_DA_rawp, zeller_polda_DA_adjustedp)
zeller_ancombc_DA_rawp <- zeller_performance$Genus[zeller_performance$ancombc_rawp < 0.05]
zeller_ancombc_DA_adjustedp <- zeller_performance$Genus[zeller_performance$ancombc_adjustedp < 0.05]
zeller_ancombc_DA_rawp <- setdiff(zeller_ancombc_DA_rawp, zeller_ancombc_DA_adjustedp)

zeller_aldex_DA_rawp <- zeller_performance$Genus[zeller_performance$aldex_rawp < 0.05]
zeller_aldex_DA_adjustedp <- zeller_performance$Genus[zeller_performance$aldex_adjustedp < 0.05]
zeller_aldex_DA_rawp <- setdiff(zeller_aldex_DA_rawp, zeller_aldex_DA_adjustedp)

zeller_result <- data.frame(POLDA_adjusted = rep('', 40),
                            POLDA_raw = rep('', 40),
                            ANCOMBC_adjusted = rep('', 40),
                            ANCOMBC_raw = rep('', 40),
                            ALDEX_adjusted = rep('', 40),
                            ALDEX_raw = rep('', 40))
zeller_result$POLDA_adjusted[1:length(zeller_polda_DA_adjustedp)] <- zeller_polda_DA_adjustedp
zeller_result$POLDA_raw[1:length(zeller_polda_DA_rawp)] <- zeller_polda_DA_rawp
zeller_result$ANCOMBC_adjusted[1:length(zeller_ancombc_DA_adjustedp)] <- zeller_ancombc_DA_adjustedp
zeller_result$ANCOMBC_raw[1:length(zeller_ancombc_DA_rawp)] <- zeller_ancombc_DA_rawp
zeller_result$ALDEX_raw[1:length(zeller_aldex_DA_rawp)] <- zeller_aldex_DA_rawp
zeller_result$ALDEX_adjusted[1:length(zeller_aldex_DA_adjustedp)] <- zeller_aldex_DA_adjustedp


library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "Zeller")
addWorksheet(wb, "Baxter")
writeData(wb, "Zeller", zeller_result, rowNames=FALSE)
writeData(wb, 'Baxter', baxter_result, rowNames=FALSE)
saveWorkbook(wb, 'real_data/metaanalysis/crc_DA_results.xlsx',
             overwrite = TRUE)

intersect(baxter_ancombc_DA_adjustedp, zeller_ancombc_DA_adjustedp)

