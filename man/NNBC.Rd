\name{NNBC}
\alias{NNBC}
\title{k nucleotide library of scheme NNBC}
\usage{
  NNBC(k)
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
NNBC(2)
}

