\name{makowski}
\alias{makowski}
\title{Diversity index according to Makowski}
\usage{
  makowski(k, libscheme)
}
\arguments{
  \item{k}{length of peptide sequences}

  \item{libscheme}{Name (character vector) or definition
  (data frame) of scheme}
}
\value{
  diversity index between 0 and 1
}
\description{
  Diversity according to Makowski is defined as ... need a
  reference to Makowski & Soares 2003 here.
}
\examples{
makowski(2, "NNN")
makowski(3, "NNK")
makowski(3, "Trimer")
}

