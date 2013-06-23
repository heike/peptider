\name{NNK}
\alias{NNK}
\title{k-nucleotide library of scheme NNK}
\usage{
  NNK(k)
}
\arguments{
  \item{k}{length of peptide sequences}
}
\value{
  list consisting of a data frame of peptide classes, size
  of class, and its probabilities, and a list of additional
  information relating to the library scheme
}
\description{
  k-nucleotide library of scheme NNK
}
\examples{
NNK(2)
head(NNK(7))
}

