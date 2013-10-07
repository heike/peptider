\name{generateCustom}
\alias{generateCustom}
\title{Generate peptide and library information for a given scheme}
\usage{
  generateCustom(scheme_name = "Custom",
    scheme_def = read.csv(file.choose()))
}
\arguments{
  \item{scheme_name}{The name of the resulting encoding
  scheme}

  \item{scheme_def}{A data frame containing encoding
  information for the scheme}
}
\value{
  TRUE upon completion of the script and output of the CSV
  files
}
\description{
  Generate peptide and library information for a given
  scheme
}
\examples{
\dontrun{
generateCustom()
generateCustom(scheme_name = "NNK_S", scheme_def = peptider:::nnk_scheme)
}
}

