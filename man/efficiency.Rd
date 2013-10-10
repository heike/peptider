\name{efficiency}
\alias{efficiency}
\title{Relative efficiency of a library}
\usage{
  efficiency(k, libscheme, N, lib = NULL)
}
\arguments{
  \item{k}{length of peptide sequences}

  \item{libscheme}{Name (character vector) or definition
  (data frame) of scheme}

  \item{N}{size of the library}

  \item{lib}{library, if null, libscheme will be used to
  create it}
}
\value{
  relative efficiency index between 0 and 1
}
\description{
  efficiency according to our paper
}
\examples{
efficiency(3, "NNN", 10^2)
efficiency(3, "NNK", 10^2)
efficiency(3, "Trimer", 10^2) ## Trimer efficiency is not 1 because of random sampling.
}

