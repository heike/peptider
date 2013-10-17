\name{scheme}
\alias{scheme}
\title{Get the specified library scheme definition}
\usage{
  scheme(name)
}
\arguments{
  \item{name}{name of the scheme as a character vector}
}
\value{
  a data frame of peptide classes, amino acids, and size of
  the classes corresponding to the selected scheme
}
\description{
  Get the specified library scheme definition
}
\examples{
scheme("NNN")
scheme("NNK")
}

