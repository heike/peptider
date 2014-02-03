\name{diversity}
\alias{diversity}
\title{Diversity according to peptides paper (Sieber)}
\usage{
diversity(k, libscheme, N, lib = NULL)
}
\arguments{
  \item{k}{length of peptide sequences}

  \item{libscheme}{Name (character vector) or definition
  (data frame) of scheme}

  \item{N}{size of the library}

  \item{lib}{library scheme}
}
\value{
Expected Diversity of the library
}
\description{
Diversity according to peptides paper (Sieber)
}
\examples{
diversity(2, "NNN", 10^3)
diversity(2, "NNK", 10^3)
}

