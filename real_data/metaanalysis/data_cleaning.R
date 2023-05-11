rm(list=ls())
library(dplyr)
library(readxl)
library(phyloseq)

gen_phylo <- function(otu, metadata){
  # otu is a matrix with the phylogeny as row names
  # metadata
  indv_IDs <- intersect(rownames(metadata), 
                        colnames(otu))
  
  otu_subset <- otu[, indv_IDs]
  metadata_subset <- metadata[indv_IDs, ]
  
  taxonomy_long <- row.names(otu_subset)
  taxonomy_break <- sapply(taxonomy_long, 
                           function(taxstring)strsplit(taxstring[1], split=';')[[1]])
  taxonomy_table <- as.matrix(taxonomy_break, row=8) |> t()
  taxonomy_table <- taxonomy_table[, seq(1, 7)] |> as.data.frame()
  colnames(taxonomy_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  rownames(otu_subset) <- sprintf("OTU_%d", seq(1, nrow(otu_subset)))
  rownames(taxonomy_table) <- sprintf("OTU_%d", 
                                      seq(1, nrow(otu_subset)))
  
  phylo_obj <- phyloseq(otu_table(as.matrix(otu_subset), taxa_are_rows=TRUE),
                             sample_data(metadata_subset),
                             tax_table(as.matrix(taxonomy_table)))
  # aggregate to genus level
  phylo_genus_obj <- tax_glom(phylo_obj, taxrank="Genus")
  # filter taxa based on prevalence
  phylo_genus_obj_filter <- filter_taxa(phylo_genus_obj, function(x) mean(x!=0)>0.01, TRUE)
  
  return(phylo_genus_obj_filter)
}


# Baxter etal from https://www.nature.com/articles/s41467-017-01973-8
crc_baxter_otu <- read.table('real_data/metaanalysis/crc_baxter_results/crc_baxter.otu_table.txt',
                             sep='\t', row.names=1, header=TRUE)
crc_baxter_metadata <- read.table('real_data/metaanalysis/crc_baxter_results/crc_baxter.metadata.txt',
                                  sep='\t', row.names=1, header=TRUE)
crc_baxter_metadata_subset <- crc_baxter_metadata[crc_baxter_metadata$DiseaseState %in% c("H", "CRC"), ]
colnames(crc_baxter_otu) <- gsub("X", "", colnames(crc_baxter_otu))
crc_baxter_phyobj <- gen_phylo(crc_baxter_otu, crc_baxter_metadata_subset)
baxter_taxa <- tax_table(crc_baxter_phyobj)
saveRDS(crc_baxter_phyobj,
        'real_data/metaanalysis/crc_baxter_results/crc_baxter_phyobj.rds')


# Zeller etal from https://www.nature.com/articles/s41467-017-01973-8
crc_zeller_otu <- read.table('real_data/metaanalysis/crc_zeller_results/crc_zeller.otu_table.txt',
                             sep='\t', row.names=1, header=TRUE)
colnames(crc_zeller_otu) <- gsub("[.]", "-", colnames(crc_zeller_otu)) 
crc_zeller_metadata <- read.table('real_data/metaanalysis/crc_zeller_results/crc_zeller.metadata.txt',
                                  sep='\t', row.names=1, header=TRUE, check.names = F)
crc_zeller_metadata_subset <- crc_zeller_metadata[crc_zeller_metadata$DiseaseState != "", ]
crc_zeller_phyobj <- gen_phylo(crc_zeller_otu, crc_zeller_metadata_subset)
zeller_taxa <- tax_table(crc_zeller_phyobj)
saveRDS(crc_zeller_phyobj,
        'real_data/metaanalysis/crc_zeller_results/crc_zeller_phyobj.rds')


# schubert etal from https://www.nature.com/articles/s41467-017-01973-8
cdi_schubert_otu <- read.table('real_data/metaanalysis/cdi_schubert_results/cdi_schubert.otu_table.txt',
                             sep='\t', row.names = 1, header=TRUE)
cdi_schubert_metadata <- read.table('real_data/metaanalysis/cdi_schubert_results/cdi_schubert.metadata.txt',
                           sep='\t', row.names = 1, header=TRUE)
cdi_schubert_phyobj <- gen_phylo(cdi_schubert_otu, cdi_schubert_metadata)

saveRDS(cdi_schubert_phyobj, 
        'real_data/metaanalysis/cdi_schubert_results/cdi_schubert_phyobj.rds')

# Vincent etal from https://www.nature.com/articles/s41467-017-01973-8
cdi_vincent_otu <- read.table('real_data/metaanalysis/cdi_vincent_v3v5_results/cdi_vincent_v3v5.otu_table.txt',
                              sep='\t', row.names = 1, header=TRUE)
cdi_vincent_metadata <- read.table('real_data/metaanalysis/cdi_vincent_v3v5_results/cdi_vincent_v3v5.metadata.txt',
                            sep='\t', row.names = 1, header=TRUE)
vincent_IDs <- gsub("-", ".", row.names(cdi_vincent_metadata))
rownames(cdi_vincent_metadata) <- vincent_IDs
cdi_vincent_phyobj <- gen_phylo(cdi_vincent_otu,
                                cdi_vincent_metadata)

cdi_vincent_phyobj_filtered <- filter_taxa(cdi_vincent_phyobj, function(x) sum(x!=0)>2, TRUE)

saveRDS(cdi_vincent_phyobj_filtered,
        'real_data/metaanalysis/cdi_vincent_v3v5_results/cdi_vincent_phyobj.rds' )

# Singh etal from https://www.nature.com/articles/s41467-017-01973-8
edd_singh_otu <- read.table('real_data/metaanalysis/edd_singh_results/edd_singh.otu_table.txt',
                            sep='\t', row.names = 1, header=TRUE)
edd_singh_metadata <- read.table('real_data/metaanalysis/edd_singh_results/edd_singh.metadata.txt',
                                 sep='\t', row.names = 1, header=TRUE)

edd_singh_phyobj <- gen_phylo(edd_singh_otu,
                              edd_singh_metadata)

saveRDS(edd_singh_phyobj,
        "real_data/metaanalysis/edd_singh_results/edd_singh_phyobj.rds")


# GEMS from https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r76
suppressMessages(library(metagenomeSeq))
library(msd16s)
data(msd16s)
gems_metadata <- phenoData(msd16s)@data
gems_taxonomy_table <- featureData(msd16s)@data
taxa_prefilter <- gems_taxonomy_table$superkingdom != 'NA'
gems_otu <- MRcounts(msd16s)
gems_otu <- gems_otu[taxa_prefilter, ]
gems_taxonomy_table <- gems_taxonomy_table[taxa_prefilter, seq(1, 7)]

gems_phyobj <- phyloseq(otu_table(as.matrix(gems_otu), taxa_are_rows = TRUE), 
                        sample_data(gems_metadata),
                        tax_table(as.matrix(gems_taxonomy_table)))

gems_genus_phyobj <- tax_glom(gems_phyobj, taxrank="genus")

gems_genus_phyobj_filter <- filter_taxa(gems_genus_phyobj, 
                                        function(x) mean(x!=0)>0.01, TRUE)

saveRDS(gems_genus_phyobj_filter,
        "real_data/metaanalysis/gems_results/gems_phyobj.rds")



# Schneider etal https://www.nature.com/articles/sdata2017152

otu_schneider <- read_excel('real_data/metaanalysis/schneider_results/Schneider_otu_table.xlsx')
taxonomy_table <- otu_schneider[, 81] 
otu_counts <- otu_schneider[, seq(2, 80)] |> as.matrix()

taxlevels <- sapply(taxonomy_table$taxonomy, 
                    function(taxstring) strsplit(taxstring, split='; ')[[1]])

level_numbers <- do.call(c, lapply(taxlevels, length))
filter_taxa <- which(level_numbers == 7)

otu_counts_subset <- otu_counts[filter_taxa, ]
taxlevels_subset <- taxlevels[filter_taxa]
taxonomy_table <- do.call(rbind, taxlevels_subset)

rownames(otu_counts_subset) <- sprintf("OTU_%d", seq(1, nrow(otu_counts_subset)))
rownames(taxonomy_table) <- sprintf("OTU_%d", seq(1, nrow(otu_counts_subset)))
colnames(taxonomy_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

metadata_schneider <- read.csv('real_data/metaanalysis/schneider_results/Schneider_metadata.csv')
metadata_schneider <- metadata_schneider[metadata_schneider$diarrhea!="", ]
rownames(metadata_schneider) <- metadata_schneider$Sample.name
metadata_schneider$Sample.name <- NULL

schneider_phyobj <- phyloseq(otu_table(otu_counts_subset, taxa_are_rows = TRUE),
                             tax_table(taxonomy_table),
                             sample_data(metadata_schneider))

schneider_genus_phyobj <- tax_glom(schneider_phyobj, taxrank="Genus")

schneider_genus_phyobj_filter <- filter_taxa(schneider_genus_phyobj, 
                                        function(x) sum(x!=0) > 3, TRUE)

saveRDS(schneider_genus_phyobj_filter,
        "real_data/metaanalysis/schneider_results/schneider_phyobj.rds")
