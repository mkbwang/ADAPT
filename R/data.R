#' Saliva samples from early childhood dental caries studies
#' 
#' A phyloseq object with 161 samples and 280 ASVs for 16S sequencing of saliva 
#' samples. The samples were collected from 
#' 12-month-old infants. 84 out of 161 children developed dental caries after 
#' 36 months old. All samples have been de-identified.
#' 
#' @usage data(ecc_saliva)
#' 
#' 
#' @format 
#' The metadata of `ecc_saliva` has two columns
#' \describe{
#'   \item{CaseStatus}{Whether the child developed dental caries after 36 months old}
#'   \item{Site}{The location of sample collection (Site 1 or Site 2).}
#' }
#' @source The original publication of this data is "Evaluating the ecological 
#' hypothesis: early life salivary microbiome assembly predicts dental caries 
#' in a longitudinal case-control study" (Blostein etal, 2022). 
#' The sequence data are available under Project number PRJNA752888.
"ecc_saliva"




#' Plaque samples from early childhood dental caries studies
#' 
#'
#' A phyloseq object with 30 samples and 610 taxa for whole genome sequencing of plaque samples. 
#' The samples were collected from children at 36, 48 or 60 months old.
#' 15 samples were from teeth with dental lesions and 15 samples were controls. 
#' The samples of cases were collected at the onset visit.
#'
#' @usage data(ecc_plaque)
#'
#' @format 
#' The metadata of `ecc_plaque` has two columns
#' \describe{
#'   \item{CaseStatus}{Whether the kids had dental caries}
#'   \item{Site}{The location of sample collection (Site 1 or Site 2).}
#' }
#' @source The original publication of this data is "Evaluating the ecological 
#' hypothesis: early life salivary microbiome assembly predicts dental caries 
#' in a longitudinal case-control study" (Blostein etal, 2022). 
#' The sequence data are available under Project number PRJNA752888.
"ecc_plaque"

