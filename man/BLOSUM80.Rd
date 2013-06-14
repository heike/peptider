\docType{data}
\name{BLOSUM80}
\alias{BLOSUM80}
\title{BLOSUM80 matrix}
\usage{
  data(BLOSUM80) Compute the number of codons for a vector
    of peptide sequences
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
  where does this matrix come from and what does it
  describe?
}
\details{
  BLOSUM80 matrix

  where does this matrix come from and what does it
  describe?
}
\examples{
codons("APE", libscheme=NNK)
codons("HENNING", libscheme=NNK)
}
\keyword{datasets}

