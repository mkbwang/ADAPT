
baxter_performance <- read.csv('real_data/metaanalysis/baxter_performance.csv',
                               row.names=1)
baxter_polda_DA_rawp <- baxter_performance$Genus[baxter_performance$polda_rawp < 0.05] |> 
  na.omit()
baxter_polda_DA_adjustedp <- baxter_performance$Genus[baxter_performance$polda_adjustedp < 0.05] |> 
  na.omit()
baxter_ancombc_DA_rawp <- baxter_performance$Genus[baxter_performance$ancombc_rawp < 0.05]
baxter_ancombc_DA_adjustedp <- baxter_performance$Genus[baxter_performance$ancombc_adjustedp < 0.05]
baxter_aldex_DA_rawp <- baxter_performance$Genus[baxter_performance$aldex_rawp < 0.05]
baxter_aldex_DA_adjustedp <- baxter_performance$Genus[baxter_performance$aldex_adjustedp < 0.05]



zeller_performance <- read.csv('real_data/metaanalysis/zeller_performance.csv',
                               row.names=1)
zeller_polda_DA_rawp <- zeller_performance$Genus[zeller_performance$polda_rawp < 0.05] |> 
  na.omit()
zeller_polda_DA_adjustedp <- zeller_performance$Genus[zeller_performance$polda_adjustedp < 0.05] |> 
  na.omit()
zeller_ancombc_DA_rawp <- zeller_performance$Genus[zeller_performance$ancombc_rawp < 0.05]
zeller_ancombc_DA_adjustedp <- zeller_performance$Genus[zeller_performance$ancombc_adjustedp < 0.05]
zeller_aldex_DA_rawp <- zeller_performance$Genus[zeller_performance$aldex_rawp < 0.05]
zeller_aldex_DA_adjustedp <- zeller_performance$Genus[zeller_performance$aldex_adjustedp < 0.05]





  