\name{ppeptide}
\alias{ppeptide}
\title{Probability of detection of a peptide sequence}
\usage{
ppeptide(x, libscheme, N)
}
\arguments{
  \item{x}{(vector) of character strings of peptide
  sequences.}

  \item{libscheme}{library scheme under which neighbors are
  being calculated.}

  \item{N}{number of valid DNA clones investigated}
}
\value{
probability of detection
}
\description{
use this function for only a few peptide sequences. Any
larger number of peptide sequences should be dealt with in
the framework of the library scheme and the detect
function.
}
\examples{
ppeptide("APE", libscheme="NNK", N=10^8)
ppeptide("HENNING", libscheme="NNK", N=10^8)
}

