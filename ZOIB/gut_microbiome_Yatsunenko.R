library(biomformat)
library(phyloseq)
library(stringr)

rm(list=ls())
folder <- "/home/wangmk/UM/Research/MDAWG/DiffRatio/RMD/Yatsunenko_data"

otu_counts <- read.table(file.path(folder, 'global_gut_otu.txt'),
                         sep="\t", header=TRUE)

metadata <- read.table(file.path(folder, 'global_gut_metadata.txt'),
                       sep='\t', header=FALSE)

taxonomy <- read.table(file.path(folder, 'global_gut_taxonomy.txt'),
                       sep='\t', header=TRUE)
