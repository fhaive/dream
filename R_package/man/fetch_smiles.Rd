% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_smiles.R
\name{fetch_smiles}
\alias{fetch_smiles}
\title{Fetch SMILES from PubChem by compound name}
\usage{
fetch_smiles(compounds)
}
\arguments{
\item{compounds}{A character vector of compound names for which SMILES strings should be fetched.}
}
\value{
A list of named SMILES strings corresponding to the input compounds.
}
\description{
This function takes a vector of compound names and retrieves their SMILES
strings from PubChem using the PubChem REST API.
}
