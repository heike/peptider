\name{codons}
\alias{codons}
\title{Compute the number of codons for a vector of peptide sequences}
\usage{
  codons(x, libscheme = NULL)
}
\arguments{
  \item{x}{(vector) of character strings of peptide
  sequences.}

  \item{libscheme}{library scheme under which neighbors are
  being calculated. this is only of importance, if
  method="dna"}
}
\value{
  vector of numbers of codons
}
\description{
  Compute the number of codons for a vector of peptide
  sequences
}
\examples{
codons("APE", libscheme=NNK)
codons("HENNING", libscheme=NNK)
}

