\name{getNofNeighbors}
\alias{getNofNeighbors}
\title{Compute the number of neighbor of degree one for a set of peptide sequences}
\usage{
  getNofNeighbors(x, blosum = 1, method = "peptide",
    libscheme = NULL)
}
\arguments{
  \item{x}{(vector) of character strings of peptide
  sequences.}

  \item{blosum}{minimal BLOSUM loading, defaults to 1 for
  positive loadings only}

  \item{method}{character string, one of "peptide" or
  "dna". This specifies the level at which the neighbors
  are calculated.}

  \item{libscheme}{library scheme under which neighbors are
  being calculated. this is only of importance, if
  method="dna"}
}
\value{
  vector of numbers of neighbors
}
\description{
  first degree neighbors - a neighbor of a peptide is
  defined as a peptide sequence that differs in at most one
  amino acid from a given sequence. Additionally, we can
  restrict neighbors to regard only those sequences that
  have a certain minimal BLOSUM loading. Use this function
  for only a few peptide sequences. Any larger number of
  peptide sequences will take too much main memory.
}
\examples{
getNofNeighbors("APE")
getNofNeighbors(c("NEAREST", "EARNEST"))
getNofNeighbors("N")
getNofNeighbors("N", method="dna", libscheme="NNK")
}

