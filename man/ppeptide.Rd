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

  \item{N}{number of DNA clones investigated}
}
\value{
  probability of detection
}
\description{
  Probability of detection of a peptide sequence
}
\examples{
ppeptide("APE", libscheme=NNK, N=10^8)
ppeptide("HENNING", libscheme=NNK, N=10^8)
}

