\name{neighbors}
\alias{neighbors}
\title{Compute number of peptide neighbors of degree 1}
\usage{
  neighbors(x, blosum = 1, type = "number",
    method = "peptide", libscheme = NULL)
}
\arguments{
  \item{x}{(vector) of character strings of peptide
  sequences.}

  \item{blosum}{minimal BLOSUM loading, defaults to 1 for
  positive loadings only}

  \item{type}{either "text" or "number", defaults to
  "number". This determines the resulting output}

  \item{method}{character string, one of "peptide" or
  "dna". This specifies the level at which the neighbors
  are calculated.}

  \item{libscheme}{library scheme under which neighbors are
  being calculated. this is only of importance, if
  method="dna"}
}
\value{
  vector of numbers of neighbors (if type="number") or list
  of sequences of the same length as x describing all
  direct neighbors of x (if type="text").
}
\description{
  first degree neighbors - a neighbor of a peptide is
  defined as a peptide sequence that differs in at most one
  amino acid from a given sequence. Additionally, we can
  restrict neighbors to regard only those sequences that
  have a certain minimal BLOSUM loading.
}
\examples{
neighbors("APE", type="text")
neighbors(c("NEAREST", "EARNEST"))
neighbors("N")
neighbors("N", method="dna", libscheme=NNK)
}

