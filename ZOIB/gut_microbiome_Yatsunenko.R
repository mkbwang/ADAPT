library(biomformat)
library(phyloseq)
library(stringr)

rm(list=ls())
folder <- "/home/wangmk/UM/Research/MDAWG/Differential_Ratio/data/study_850"

modify_string <- function(mystring){
  if (!grepl('.', mystring, fixed=TRUE)){ # append .1 at the end
    newstring <- paste(mystring, '1', sep=".")
  } else{
    newstring <- mystring
  }
  if(substr(newstring, 1, 1) == 'k'){
    newstring <- sub('^.', 'h', newstring)
  }
  return(newstring)
}

## try to align the individual names in the otumat column with the sample metadata sheet
otu_table <- read_biom(file.path(folder, 'processed_data', '421_otu_table.biom'))
otumat <- as(biom_data(otu_table), "matrix")
indv_names <- colnames(otumat) %>% unname()
indv_names <- str_sub(indv_names, 5)
indv_names <- sapply(indv_names, modify_string) %>% unname()
colnames(otumat) <- indv_names


sample_info <- read.csv(file.path(folder, 'supplementary_table', 'individuals.csv'))
sample_info$Sample_Identifier[sample_info$Sample_Identifier == 'h10A.2'] <- 'h10A.1'
sample_info$Sample_Identifier[sample_info$Sample_Identifier == 'h257M.2'] <- 'h257M.1'
sample_info$Sample_Identifier[sample_info$Sample_Identifier == 'h37A.2'] <- 'h37A.1'
sample_info$Sample_Identifier[sample_info$Sample_Identifier == 'h37B.2'] <- 'h37B.1'
row.names(sample_info) <- sapply(sample_info$Sample_Identifier, modify_string)
sample_info$Sample_Identifier <- NULL

sample_info <- sample_info[colnames(otumat), ]

# OTU metadata is available in the biom file
OTU_metadata <- observation_metadata(otu_table)


# set up phyloseq object
OTU_object <- otu_table(otumat, taxa_are_rows = TRUE)
TAX_object <- tax_table(OTU_metadata)
taxa_names(TAX_object) <- row.names(otumat)
SAMPLE_object <- sample_data(sample_info)

phyloseq_object <- phyloseq(OTU_object, TAX_object, SAMPLE_object)

saveRDS(phyloseq_object, file=file.path(folder, 'Yatsunenko_phyloseq.rds'))
