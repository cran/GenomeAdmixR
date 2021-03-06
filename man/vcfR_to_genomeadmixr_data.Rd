% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_input_data.R
\name{vcfR_to_genomeadmixr_data}
\alias{vcfR_to_genomeadmixr_data}
\title{function to convert a vcfR object to genome_admixr_data}
\usage{
vcfR_to_genomeadmixr_data(
  vcfr_object,
  chosen_chromosome,
  number_of_snps = NA,
  random_snps = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{vcfr_object}{result of vcfR::read.vcfR}

\item{chosen_chromosome}{chromosome of choice}

\item{number_of_snps}{number of snps to be loaded from the vcf file, default
is to load all snps}

\item{random_snps}{if a subset of all snps has to be taken, should these
be sampled sequentially (e.g. the first 100 snps) or randomly (100 randomly
sampled snps) (examples are for 'number_of_snps' = 100).}

\item{verbose}{if true, print progress bar}
}
\value{
genomeadmixr_data object ready for simulate_admixture_data
}
\description{
function to convert a vcfR object to genome_admixr_data
}
