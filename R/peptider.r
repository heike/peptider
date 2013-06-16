
#' Nucleotide library of scheme NNN
#'
#' @param k length of peptide sequences
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, 
#' and a list of additional information relating to the library scheme
#' @export
#' @examples
#' NNN(2)
#' head(NNN(7))
NNN <- function(k) {
#   seq <- make.RV(c("A","B","C", "D", "E"), c(3*6, 5*4, 1*3, 8*2, 2*1))
#   d <- make.RV(c("A","B","C", "D", "E"), c(3, 5, 1, 8, 2))
#   
#   d7 <- multN(d,k)
#   seq7 <- multN(seq,k)
#   
#   di <- round(probs(d7)*19^k,0)
#   pi <- probs(seq7)
#   
#   list(data = data.frame(class = as.vector(d7), di = di, probs = pi), 
#        info =list(nucleotides = 64, valid = 59, 
#                   scheme = data.frame(class=c("A", "B", "C", "D", "E", "Z"), 
#                                   aacids=c("SLR", "AGPTV", "I", "DEFHKNQY", "MW", "C*"), 
#                                   s=c(3,5,1,8,2,2), c=c(6,4,3,2,1,2.5))))
  libBuild(k, libscheme=nnn_scheme)
}

#' Nucleotide library  scheme NNN
#'
#' @name nnn_scheme
#' @title Nucleotide library scheme NNN
#' @description This data set contains descriptions of amino acid classes under the NNN library scheme.
#' @docType data
#' @usage libBuild(1, libscheme="nnn_scheme")

nnn_scheme <- data.frame(class=c("A", "B", "C", "D", "E", "Z"), 
                            aacids=c("SLR", "AGPTV", "I", "DEFHKNQY", "MW", "C*"), 
                            s=c(3,5,1,8,2,2), c=c(6,4,3,2,1,2.5), stops=3)


#' Nucleotide library of scheme Trimer
#'
#' @param k length of peptide sequences
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, 
#' and a list of additional information relating to the library scheme
#' @export
#' @examples
#' Trimer(2)
#' head(Trimer(7))
Trimer <- function(k) {
#   seq <- make.RV(c("A"), c(19))
#   d <- make.RV(c("A"), c(19))
#   
#   d7 <- multN(d,k)
#   seq7 <- multN(seq,k)
#   
#   di <- round(probs(d7)*19^k,0)
#   pi <- probs(seq7)
#   
#   list(data = data.frame(class = as.vector(d7), di = di, probs = pi), 
#        info =list(nucleotides = 19, valid = 19, 
#                   scheme = data.frame(class=c("A"), 
#                                   aacids=c("SLRAGPTVIDEFHKNQYMW"), 
#                                   s=c(19), c=c(1))))
  libBuild(k, libscheme=trimer_scheme)
}

#' Nucleotide library  scheme Trimer
#'
#' Only valid peptide sequences are created. Each sequence appears with the same probability.
#' @name trimer_scheme
#' @title Nucleotide library scheme Trimer
#' @description This data set contains descriptions of amino acid classes under the Trimer library scheme.
#' @docType data
#' @usage libBuild(1, libscheme="trimer_scheme")

trimer_scheme <- data.frame(class=c("A", "Z"), 
                            aacids=c("SLRAGPTVIDEFHKNQYMW", "*"), 
                            s=c(19, 0), c=c(1, 0), stops=0)






#' Diversity index according to Makowski
#'
#' Diversity according to Makowski is defined as ... need a reference to Makowski & Soares 2003 here. 
#' @param k length of peptide sequences
#' @param libscheme function 
#' @return diversity index between 0 and 1
#' @export
#' @examples
#' makowski(2, NNN)
#' makowski(3, NNK)
#' makowski(3, Trimer)
makowski <- function(k, libscheme) {
  dframe <- libscheme(k)$data
  with(dframe, 1/(19^k*sum(probs^2/di)))
}

#' Coverage as expected number of peptides given all possible peptides
#'
#' Coverage of library of size N given random sampling from the pool of all possible peptides according 
#' to probabilities determined according to the library scheme.
#' @param k length of peptide sequences
#' @param libscheme function 
#' @param N size of the library 
#' @param lib library, if null, libscheme will be used to create it
#' @return coverage index between 0 and 1
#' @export
#' @examples
#' coverage(2, NNN, 10^3)
#' coverage(2, NNK, 10^3)
#' coverage(2, Trimer, 10^3) ## Trimer coverage is not 1 because of random sampling.
coverage <- function(k, libscheme, N, lib=NULL) {
  if (is.null(lib)) lib <- libscheme(k)
  libdata <- lib$data
  
  initialloss <- (1-(lib$info$valid/lib$info$nucleotides)^k)
  libdata$expected <- libdata$probs*N*(1-initialloss)
  libdata$z <- with(libdata, di*(1-exp(-expected/di)))
  
  with(libdata, min(sum(z)/19^k,1))
}

#' Relative efficiency of a library
#'
#' efficiency according to our paper
#' @param k length of peptide sequences
#' @param libscheme function 
#' @param N size of the library 
#' @param lib library, if null, libscheme will be used to create it
#' @return relative efficiency index between 0 and 1
#' @export
#' @examples
#' efficiency(3, NNN, 10^2)
#' efficiency(3, NNK, 10^2)
#' efficiency(3, Trimer, 10^2) ## Trimer efficiency is not 1 because of random sampling.
efficiency <- function(k, libscheme, N, lib=NULL) {
  if (is.null(lib)) lib <- libscheme(k)
  libdata <- lib$data
  
  initialloss <- (1-(lib$info$valid/lib$info$nucleotides)^k)
  libdata$expected <- libdata$probs*N*(1-initialloss)
  libdata$z <- with(libdata, di*(1-exp(-expected/di)))
  
  with(libdata, sum(z)/min(19^k,N))
}


#' Nucleotide library of scheme NNK
#'
#' @param k length of peptide sequences
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, 
#' and a list of additional information relating to the library scheme
#' @export
#' @examples
#' NNK(2)
#' head(NNK(7))
NNK <- function(k) {
#   seq <- make.RV(c("A","B","C"), c(3*3,5*2,11))
#   d <- make.RV(c("A","B","C"), c(3,5,11))
#   
#   d7 <- multN(d,k)
#   seq7 <- multN(seq,k)
#   
#   di <- round(probs(d7)*19^k,0)
#   pi <- probs(seq7)
#   
#   list(data=data.frame(class = as.vector(d7), di = di, probs = pi),
#        info=list(nucleotides=32, valid=30,
#                  scheme=data.frame(class=c("A", "B", "C", "Z"),
#                               aacid=c("SLR", "AGPTV", "DEFHIKMNQWY", "C*"),
#                                s=c(3,5,11,2),
#                                c=c(3,2,1,1))))
  libBuild(k, nnk_scheme)
}

#' Nucleotide library  scheme NNK
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, G, or T.
#' @name nnk_scheme
#' @title Nucleotide library scheme NNK
#' @description This data set contains descriptions of amino acid classes under the NNK library scheme.
#' @docType data
#' @usage libBuild(1, libscheme=nnk_scheme)

nnk_scheme <- data.frame(class=c("A", "B", "C", "Z"),
                         aacid=c("SLR", "AGPTV", "DEFHIKMNQWY", "C*"),
                         s=c(3,5,11,2),
                         c=c(3,2,1,1), stops=2)


#' Nucleotide library of scheme NNS
#'
#' @param k length of peptide sequences
#' @return data frame of peptide classes, size of class, and its probability
#' @export
#' @examples
#' NNS(2)
#' head(NNS(7))
NNS <- function(k) {
  ## from a codon to amino acid encoding point of view NNK and NNS are identical.
  NNK(k)
}

#' Nucleotide library of scheme NNB
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, G, or T.
#' @param k length of peptide sequences
#' @return data frame of peptide classes, size of class, and its probability
#' @export
#' @examples
#' NNB(2)
#' head(NNB(7))
NNB <- function(k) {
#   seq <- make.RV(c("A","B","C", "D", "E"), c(1*5, 2*4, 5*3, 6*2, 5*1))
#   d <- make.RV(c("A","B","C", "D", "E"), c(1,2,5,6,5))
#   
#   d7 <- multN(d,k)
#   seq7 <- multN(seq,k)
#   
#   di <- round(probs(d7)*19^k,0)
#   pi <- probs(seq7)
#   
#   list(data=data.frame(class = as.vector(d7), di = di, probs = pi),
#        info=list(nucleotides=48, valid=45,
#                  scheme=data.frame(class=c("A", "B", "C", "D", "E", "Z"),
#                                aacid=c("S", "LR", "AGPTV", "DFHINY", "EKMQW", "C*"),
#                                s=c(1,2,5,6,5,2),
#                                c=c(5,4,3,2,1, 1.5))))
  libBuild(k, libscheme=nnb_scheme)
}

#' Nucleotide library  scheme NNB
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, G, or T.
#' @name nnb_scheme
#' @title Nucleotide library scheme NNB
#' @description This data set contains descriptions of amino acid classes under the NNB library scheme.
#' @docType data
#' @usage libBuild(1, libscheme="nnb_scheme")

nnb_scheme <- data.frame(class=c("A", "B", "C", "D", "E", "Z"),
              aacid=c("S", "LR", "AGPTV", "DFHINY", "EKMQW", "*C"),
              s=c(1,2,5,6,5,2),
              c=c(5,4,3,2,1, 1.5), stops=3)

#' Build peptide library of k-length sequences according to specified scheme
#' 
#' some more explanation here
#' @param k length of peptide sequences
#' @param libscheme library scheme specifying classes of amino acids according to number of encodings
#' last class is reserved for stop tags and other amino acids we are not interested in. 
#' @return library and library scheme used
#' @export
#' @examples
#' lib <- libBuild(4, libscheme=nnb_scheme)
 
libBuild <- function(k, libscheme) {
  libscheme$class <- as.character(libscheme$class)
  seq <- with(libscheme[-nrow(libscheme),], make.RV(class, s*c))
  d <- with(libscheme[-nrow(libscheme),], make.RV(class, s))
  
  d7 <- multN(d,k)
  seq7 <- multN(seq,k)
  di <- with(libscheme, round(probs(d7)*sum(s[-length(unique(class))])^k,0))
  pi <- probs(seq7)
  list(data=data.frame(class = as.vector(d7), di = di, probs = pi),
       info=list(nucleotides=sum(with(libscheme, s*c)), 
                 valid=with(libscheme, sum(s*c)-stops[1]),
                 scheme=libscheme))
}


#' Detection probability in a single library of size N
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, G, or T.
#' @param lib library used in experiment, defaults to NNK(7)
#' @param size size of the library, defaults to 10^8
#' @return vector of detection probabilities for peptide sequences in each class
#' @export
#' @examples
#' summary(detect())
#'
#' require(ggplot2)
#' lib = NNK(7)
#' qplot(detect(lib), weight=di, geom="histogram")

detect <- function(lib = NNK(7), size = 10^8) {
  with(lib$data, 1 - exp(-size*probs/di))
}


getNeighborOne <- function(x, blosum=1) {
  data(BLOSUM80)
  replacements <- llply(strsplit(x,""), function(y) {
    llply(y, function(z) {
      as.character(subset(BLOSUM80, (AA1 == z) & (Blosum >= blosum)& (AA2 != z))$AA2 )
    })
  })[[1]]
  neighbors <- NULL
    for (i in 1:nchar(x)) {
      neighbors <- c(neighbors, paste(substr(x, 1,i-1), replacements[[i]], substr(x, i+1, nchar(x)), sep=""))
    }
  # check that all neighbors have the correct length
  idx <- which(nchar(neighbors) != nchar(x))
  if (length(idx)>0) neighbors <- neighbors[-idx]
  neighbors <- unique(c(x, neighbors))
  
  return(neighbors)
}

#' Find all neighbors of degree one for a set of peptide sequences
#' 
#' first degree neighbors - a neighbor of a peptide is defined as a peptide sequence that differs in at most one amino acid from a given sequence. 
#' Additionally, we can restrict neighbors to regard only those sequences that have a certain minimal BLOSUM loading. 
#' @param x (vector) of character strings of  peptide sequences.
#' @param blosum minimal BLOSUM loading, defaults to 1 for positive loadings only
#' @return list of neighbor sequences
#' @export
#' @examples
#' getNeighbors("APE")
#' getNeighbors(c("HI", "APE"))
#' getNeighbors(c("HI", "EARNEST", "APE"), blosum=3)
#' ## degree 2 neighbors:
#' unique(unlist(getNeighbors(getNeighbors("APE"))))
getNeighbors <- function(x, blosum=1) {
  if (length(x) == 1) return(getNeighborOne(x, blosum))
  llply(x, getNeighborOne)
}

getNofNeighborsOne <- function(x, blosum = 1, method="peptide", libscheme=NULL) {
  data(BLOSUM80)
  replacements <- llply(strsplit(x,""), function(y) {
    llply(y, function(z) {
      as.character(subset(BLOSUM80, (AA1 == z) & (Blosum >= blosum) & (AA2 != z))$AA2)
    })
  })[[1]]
  if (method == "peptide") return(length(unlist(replacements))+1)
  
  stopifnot(!is.null(libscheme))
  lib <- libscheme(1)$info$scheme
   
  dnas <- sum(unlist(llply(unlist(replacements), function(w) { 
    lib[grep(w, lib$aacid),"c"]
  })))
  dnas <- dnas + sum(unlist(llply(unlist(strsplit(x, split="")), function(w) { 
    lib[grep(w, lib$aacid),"c"]
  })))
  return(dnas)
}

#' Compute the number of neighbor of degree one for a set of peptide sequences
#' 
#' first degree neighbors - a neighbor of a peptide is defined as a peptide sequence that differs in at most one amino acid from a given sequence. 
#' Additionally, we can restrict neighbors to regard only those sequences that have a certain minimal BLOSUM loading. 
#' @param x (vector) of character strings of  peptide sequences.
#' @param blosum minimal BLOSUM loading, defaults to 1 for positive loadings only
#' @param method character string, one of "peptide" or "dna". This specifies the level at which the neighbors are calculated.
#' @param libscheme library scheme under which neighbors are being calculated. this is only of importance, if method="dna"
#' @return vector of numbers of neighbors 
#' @export
#' @examples
#' getNofNeighbors("APE")
#' getNofNeighbors(c("NEAREST", "EARNEST"))
#' getNofNeighbors("N")
#' getNofNeighbors("N", method="dna", libscheme=NNK)
getNofNeighbors <- function(x, blosum = 1, method="peptide", libscheme=NULL) {
  data(BLOSUM80)
  if (length(x) == 1) return(getNofNeighborsOne(x, blosum, method, libscheme))

  return(llply(x, getNofNeighborsOne, blosum, method, libscheme))
}




#' Compute the number of codons for a vector of peptide sequences
#' 
#' @param x (vector) of character strings of  peptide sequences.
#' @param libscheme library scheme under which neighbors are being calculated. this is only of importance, if method="dna"
#' @return vector of numbers of codons 
#' @export
#' @examples
#' codons("APE", libscheme=NNK)
#' codons("HENNING", libscheme=NNK)
codons <- function(x, libscheme=NULL) {
  if (length(x) == 1) return(codonsOne(x, libscheme))

  unlist(llply(x, codonsOne, libscheme=libscheme))
}

codonsOne <- function(x, libscheme) {
  stopifnot(!is.null(libscheme))
  lib <- libscheme(1)$info$scheme

  prod(unlist(llply(strsplit(x, split="")[[1]], function(w) { 
    lib[grep(w, lib$aacid),"c"]
  })))
}

#' Probability of detection of a peptide sequence 
#' 
#' @param x (vector) of character strings of  peptide sequences.
#' @param libscheme library scheme under which neighbors are being calculated. 
#' @param N number of valid DNA clones investigated
#' @return probability of detection
#' @export
#' @examples
#' ppeptide("APE", libscheme=NNK, N=10^8)
#' ppeptide("HENNING", libscheme=NNK, N=10^8)

ppeptide <- function(x, libscheme, N) {
  n <- sum(codons(x, libscheme=libscheme))
  Max <- libscheme(1)$info$valid^nchar(as.character(x[1]))
  1 - exp(-N*n/Max)
}


#' BLOSUM80 matrix
#' 
#' where does this matrix come from and what does it describe?
#' @name BLOSUM80
#' @title BLOSUM80 matrix
#' @description where does this matrix come from and what does it describe?
#' @docType data
#' @usage data(BLOSUM80)
data(BLOSUM80)