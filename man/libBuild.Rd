\name{libBuild}
\alias{libBuild}
\title{Build peptide library of k-length sequences according to specified scheme}
\usage{
  libBuild(k, libscheme)
}
\arguments{
  \item{k}{length of peptide sequences}

  \item{libscheme}{library scheme specifying classes of
  amino acids according to number of encodings last class
  is reserved for stop tags and other amino acids we are
  not interested in.}
}
\value{
  library and library scheme used
}
\description{
  some more explanation here
}
\examples{
user_scheme <- data.frame(class=c("A", "B", "C", "Z"),
                          aacid=c("SLR", "AGPTV", "CDEFHIKMNQWY", "*"),
                          c=c(3,2,1,1))
user_library <- libBuild(3, user_scheme)
}

