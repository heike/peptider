\name{detect}
\alias{detect}
\title{Detection probability in a single library of size N}
\usage{
  detect(lib = NNK(7), size = 10^8)
}
\arguments{
  \item{lib}{library used in experiment, defaults to
  NNK(7)}

  \item{size}{size of the library, defaults to 10^8}
}
\value{
  vector of detection probabilities for peptide sequences
  in each class
}
\description{
  The last DNA nucleus in the sequence is restricted to be
  one of C, G, or T.
}
\examples{
summary(detect())

require(ggplot2)
lib = NNK(7)
qplot(detect(lib, size=10^8), weight=di, geom="histogram", data=lib$data)
}

