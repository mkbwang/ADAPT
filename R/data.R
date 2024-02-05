#' 16S sequencing of saliva samples collected from 12-month-old infants
#'
#' A phyloseq object with 161 samples and 280 ASVs. The samples are collected from 12-month-old infants' saliva.
#' 84 out of 161 children developed dental caries after 36 months old. The children were from Pennsylvania and West Virginia.
#'
#' @format ## `ecc16s_12`
#' The metadata has three columns
#' \describe{
#'   \item{CaseEver}{Whether the child developed dental caries after 36 months old}
#'   \item{geo_loc_name}{Where the child lives (Pennsylvania or West Virginia)}
#' }
#' @source The original publication of this data is "Evaluating the ecological hypothesis: early life salivary microbiome assembly predicts dental caries in a longitudinal case-control study" (Blostein etal, 2022). The sequence data are available under Project number PRJNA752888.
"ecc16s_12"






#' Whole genome sequencing of plaque samples collected from infants after 36 months old
#'
#' A phyloseq object with 30 samples and 610 taxa. The samples are collected from plaque of children at 36, 48 or 60 months old
#' 15 samples were from teeth with dental lesions and 15 samples were controls.
#'
#' @format ## `ecc_metag`
#' The metadata has three columns
#' \describe{
#'   \item{CaseStatus}{Whether the kids had dental caries}
#'   \item{Site}{Where the child lives (Pennsylvania or West Virginia)}
#' }
#' @source The original publication of this data is "Evaluating the ecological hypothesis: early life salivary microbiome assembly predicts dental caries in a longitudinal case-control study" (Blostein etal, 2022). The sequence data are available under Project number PRJNA752888.
"ecc_metag"

