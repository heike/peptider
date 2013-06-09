\name{NNS}
\alias{NNS}
\title{Nucleotide library of scheme NNS}
\usage{
  NNS(k)
}
\arguments{
  \item{k}{length of peptide sequences}
}
\value{
  data frame of peptide classes, size of class, and its
  probability
}
\description{
  Nucleotide library of scheme NNS
}
\examples{
NNS(2)
head(NNS(7))
}

