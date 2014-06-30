\name{calcNeighbors}
\alias{calcNeighbors}
\title{Calculate neighborhood distribution}
\usage{
calcNeighbors(lib, k)
}
\arguments{
  \item{lib}{library scheme}

  \item{k}{length of the peptide sequences}
}
\value{
dataset of peptide sequences: AA are amino acid sequences,
c0 are codons for self representation, cr is the ratio of
#neighbors in first degree neighborhood (not counting self
representations) and #codons in self representation N1 is
the number of neighbors in codon representation (including
self representation)
}
\description{
Calculate distribution of neighbors under library scheme
lib for peptide sequences of length k.
}
\examples{
calcNeighbors("NNK", 2)
calcNeighbors("Trimer", 2)
}

