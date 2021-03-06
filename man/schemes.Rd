% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptider.r
\docType{data}
\name{schemes}
\alias{schemes}
\title{Built-in library schemes for peptider}
\usage{
data(schemes)
}
\description{
This data set contains descriptions of amino acid classes several commonly used library schemes: NNN, NNB, NNK, 20/20, and variations of each in which Cysteine is not considered a viable amino acid.
}
\details{
Built-in library schemes


The schemes are defined as:

NNN: All four bases (\"N\" = G/A/T/C) possible at all three positions in the codon.
NNB: All four bases in the first two codon positions possible, the third position is restricted to G, T or C (= \"B\")
NNK/S: All four bases in the first two codon positions possible, the third position is restricted to G/T (= \"K\") or two C/G (= \"S\").
2020: 20/20 describes the concept that DNA is assembled from prefabricated trimeric building blocks. This allows the generation of libraries from a predefined set of codons and thereby complete exclusion of Stop codons and other unwanted codons.
NNN (-C): NNN with Cysteine ignored.
NNB (-C): NNB with Cysteine ignored.
NNK/SC (-C): NNK/S with Cysteine ignored.
2020 (-C): 20/20 with Cysteine ignored.

The schemes differ in the number of used codons, ranging from 64 (NNN), 48 (NNB), 32 (NNK/S) to 20 or less (20/20). Coding schemes that allow varying ratios of codons/amino acid, result in libraries biased towards amino acids which are encoded more often. Further, the number of Stop codons that can lead to premature termination of the peptide sequence influences the performance of the library.
}

