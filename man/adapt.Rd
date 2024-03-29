% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ADAPT.R
\name{adapt}
\alias{adapt}
\title{Pooling Tobit Models for microbiome differential abundance analysis}
\usage{
adapt(
  otu_table,
  metadata,
  covar,
  adjust = NULL,
  prevalence_cutoff = 0.05,
  depth_cutoff = 1000,
  taxa_are_rows = TRUE,
  boot = TRUE,
  boot_replicate = 2000,
  n_boot_taxa = 500,
  alpha = 0.05
)
}
\arguments{
\item{otu_table}{microbiome abundance table generated from 16S rRNA sequencing or shotgun metagenomic sequencing. All entries are integers}

\item{metadata}{sample metadata dataframe}

\item{covar}{the name of the covariate of interest}

\item{adjust}{the names of confounders to adjust for, NULL if no adjustment}

\item{prevalence_cutoff}{taxa whose prevalence are smaller than the cutoff will be excluded from analysis}

\item{depth_cutoff}{a sample would be discarded if its sequencing depth is smaller than the threshold}

\item{taxa_are_rows}{whether the microbiome taxa are on the column or the row of the abundance table}

\item{boot}{whether to use bootstrap to estimate scaling factors for the test statistics, default TRUE}

\item{boot_replicate}{number of bootstrap replicates for each taxon to estimate scaling factors}

\item{n_boot_taxa}{number of taxa to apply bootstrap}

\item{alpha}{the cutoff of the adjusted p values}
}
\value{
the reference taxa set, the identified differentially abundant taxa and the p values for all the taxa
}
\description{
Pooling Tobit Models for microbiome differential abundance analysis
}
