\name{NNN}
\alias{NNN}
\title{k-nucleotide library of scheme NNN}
\usage{
  NNN(k)
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
  k-nucleotide library of scheme NNN
}
\examples{
NNN(2)
head(NNN(7))
}

