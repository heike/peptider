\name{Trimer}
\alias{Trimer}
\title{Nucleotide library of scheme Trimer}
\usage{
  Trimer(k)
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
  Nucleotide library of scheme Trimer
}
\examples{
Trimer(2)
head(Trimer(7))
}

