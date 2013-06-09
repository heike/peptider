\name{NNB}
\alias{NNB}
\title{Nucleotide library of scheme NNB}
\usage{
  NNB(k)
}
\arguments{
  \item{k}{length of peptide sequences}
}
\value{
  data frame of peptide classes, size of class, and its
  probability
}
\description{
  The last DNA nucleus in the sequence is restricted to be
  one of C, G, or T.
}
\examples{
NNB(2)
head(NNB(7))
}

