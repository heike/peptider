\name{generateCustomProbs}
\alias{generateCustomProbs}
\title{For a given scheme, generate a dataset with the peptide probabilities}
\usage{
generateCustomProbs(scheme_def, k = 6:10, n = 6:14)
}
\arguments{
  \item{scheme_def}{definition of the custom scheme}

  \item{k}{peptide lengths to include}

  \item{n}{exponents of the library size to include}
}
\value{
A data frame of peptide probabilities
}
\description{
For a given scheme, generate a dataset with the peptide
probabilities
}

