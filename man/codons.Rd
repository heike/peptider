\name{codons}
\alias{codons}
\title{Compute the number of codons for a vector of peptide sequences}
\usage{
  codons(x, libscheme = NULL, flag = FALSE)
}
\arguments{
  \item{x}{(vector) of character strings of peptide
  sequences.}

  \item{libscheme}{library scheme under which neighbors are
  being calculated. this is only of importance, if
  method="dna"}

  \item{flag}{internal use only: Set to true if calling
  this from another function}
}
\value{
  vector of numbers of codons
}
\description{
  use this function for only a few peptide sequences. Any
  larger number of peptide sequences should be dealt with
  in the framework of the library scheme and the detect
  function.
}
\examples{
codons("APE", libscheme="NNK")
codons("HENNING", libscheme="NNK")
}

