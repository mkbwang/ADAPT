library(phyloseq)
library(mia)
library(foreach)
library(doParallel)
library(dplyr)
cores=parallelly::availableCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)



