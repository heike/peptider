
#' k-nucleotide library of scheme NNN
#'
#' @param k length of peptide sequences
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, 
#' and a list of additional information relating to the library scheme
#' @export
#' @examples
#' NNN(2)
NNN <- function(k) {

  libBuild(k, libscheme=nnn_scheme)
}

#' k-nucleotide library of scheme NNNC
#'
#' @param k length of peptide sequences
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, 
#' and a list of additional information relating to the library scheme
#' @export
#' @examples
#' NNNC(2)
NNNC <- function(k) {
    
    libBuild(k, libscheme=nnnc_scheme)
}

#' Nucleotide library  scheme NNN
#'
#' @name nnn_scheme
#' @title Nucleotide library scheme NNN
#' @description This data set contains descriptions of amino acid classes under the NNN library scheme.
#' @docType data
#' @usage libBuild(1, libscheme=nnn_scheme)

nnn_scheme <- data.frame(class=c("A", "B", "C", "D", "E", "Z"), 
                            aacids=c("SLR", "AGPTV", "I", "DEFHKNQY", "MW", "C*"), 
                            c=c(6,4,3,2,1,2.5))

#' Nucleotide library scheme NNNC
#'
#' @name nnnc_scheme
#' @title Nucleotide library scheme NNNC
#' @description This data set contains descriptions of amino acid classes under the NNN library scheme with C.
#' @docType data
#' @usage libBuild(1, libscheme=nnnc_scheme)

nnnc_scheme <- data.frame(class=c("A", "B", "C", "D", "E", "Z"), 
                         aacids=c("SLR", "AGPTV", "I", "DEFHKNQYC", "MW", "*"), 
                         c=c(6,4,3,2,1,3))


#' k-nucleotide library of scheme Trimer
#'
#' @param k length of peptide sequences
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, 
#' and a list of additional information relating to the library scheme
#' @export
#' @examples
#' Trimer(2)
Trimer <- function(k) {

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
                            c=c(1, 0))






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
  info <- libscheme(1)$info$scheme
  numAA <- sum(info$s[-nrow(info)]) 
  with(dframe, 1/(numAA^k*sum(probs^2/di)))
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
  
  s_count <- sum(subset(lib$info$scheme, class != "Z")$s)
  
  with(libdata, min(sum(z)/s_count^k,1))
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
  
  s_count <- sum(subset(lib$info$scheme, class != "Z")$s)
  
  with(libdata, min(s_count^k,sum(z))/N)
}


#' k-nucleotide library of scheme NNK
#'
#' @param k length of peptide sequences
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, 
#' and a list of additional information relating to the library scheme
#' @export
#' @examples
#' NNK(2)
NNK <- function(k) {

  libBuild(k, nnk_scheme)
}

#' k-nucleotide library of scheme NNKC
#'
#' @param k length of peptide sequences
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, 
#' and a list of additional information relating to the library scheme
#' @export
#' @examples
#' NNKC(2)
NNKC <- function(k) {
    
    libBuild(k, nnkc_scheme)
}

#' Nucleotide library  scheme NNK
#'
#' The last DNA nucleus in the sequence is restricted to be one of G, or T.
#' @name nnk_scheme
#' @title Nucleotide library scheme NNK
#' @description This data set contains descriptions of amino acid classes under the NNK library scheme.
#' @docType data
#' @usage libBuild(1, libscheme=nnk_scheme)

nnk_scheme <- data.frame(class=c("A", "B", "C", "Z"),
                         aacid=c("SLR", "AGPTV", "DEFHIKMNQWY", "C*"),
                         c=c(3,2,1,1))

#' Nucleotide library  scheme NNKC
#'
#' The last DNA nucleus in the sequence is restricted to be one of G, or T.
#' @name nnkc_scheme
#' @title Nucleotide library scheme NNKC
#' @description This data set contains descriptions of amino acid classes under the NNKC library scheme.
#' @docType data
#' @usage libBuild(1, libscheme=nnkc_scheme)

nnkc_scheme <- data.frame(class=c("A", "B", "C", "Z"),
                         aacid=c("SLR", "AGPTV", "DEFHIKMNQWYC", "*"),
                         c=c(3,2,1,1))

#' Nucleotide library  scheme NNS
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, or G.
#' @name nns_scheme
#' @title Nucleotide library scheme NNS
#' @description This data set contains descriptions of amino acid classes under the NNS library scheme.
#' @docType data
#' @usage libBuild(1, libscheme=nns_scheme)
nns_scheme <- nnk_scheme

#' Nucleotide library  scheme NNSC
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, or G.
#' @name nnsc_scheme
#' @title Nucleotide library scheme NNSC
#' @description This data set contains descriptions of amino acid classes under the NNSC library scheme.
#' @docType data
#' @usage libBuild(1, libscheme=nnsc_scheme)
nnsc_scheme <- nnkc_scheme

#' k nucleotide library of scheme NNS
#'
#' @param k length of peptide sequences
#' @return data frame of peptide classes, size of class, and its probability
#' @export
#' @examples
#' NNS(2)
NNS <- function(k) {
  ## from a codon to amino acid encoding point of view NNK and NNS are identical.
  NNK(k)
}

#' k nucleotide library of scheme NNSC
#'
#' @param k length of peptide sequences
#' @return data frame of peptide classes, size of class, and its probability
#' @export
#' @examples
#' NNSC(2)
NNSC <- function(k) {
    ## from a codon to amino acid encoding point of view NNKC and NNSC are identical.
    NNKC(k)
}

#' k nucleotide library of scheme NNB
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, G, or T.
#' @param k length of peptide sequences
#' @return data frame of peptide classes, size of class, and its probability
#' @export
#' @examples
#' NNB(2)
NNB <- function(k) {

  libBuild(k, libscheme=nnb_scheme)
}

#' k nucleotide library of scheme NNBC
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, G, or T.
#' @param k length of peptide sequences
#' @return data frame of peptide classes, size of class, and its probability
#' @export
#' @examples
#' NNBC(2)
NNBC <- function(k) {
    
    libBuild(k, libscheme=nnbc_scheme)
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
              c=c(5,4,3,2,1, 1.5))

#' Nucleotide library  scheme NNBC
#'
#' The last DNA nucleus in the sequence is restricted to be one of C, G, or T.
#' @name nnbc_scheme
#' @title Nucleotide library scheme NNBC
#' @description This data set contains descriptions of amino acid classes under the NNBC library scheme with C.
#' @docType data
#' @usage libBuild(1, libscheme="nnbc_scheme")

nnbc_scheme <- data.frame(class=c("A", "B", "C", "D", "E", "Z"),
                         aacid=c("S", "LR", "AGPTV", "DFHINYC", "EKMQW", "*"),
                         c=c(5,4,3,2,1,1))

#' Build peptide library of k-length sequences according to specified scheme
#' 
#' some more explanation here
#' @param k length of peptide sequences
#' @param libscheme library scheme specifying classes of amino acids according to number of encodings
#' last class is reserved for stop tags and other amino acids we are not interested in. 
#' @return library and library scheme used
#' @examples
#' user_scheme <- data.frame(class=c("A", "B", "C", "Z"),
#'                           aacid=c("SLR", "AGPTV", "CDEFHIKMNQWY", "*"),
#'                           c=c(3,2,1,1))
#' user_library <- libBuild(3, user_scheme)                        
#' @export
 
libBuild <- function(k, libscheme) {
  libscheme$class <- as.character(libscheme$class)
  libscheme$s <- nchar(as.character(libscheme$aacid))
  seq <- with(libscheme[-nrow(libscheme),], make.RV(class, s*c))
  d <- with(libscheme[-nrow(libscheme),], make.RV(class, s))
  
  d7 <- multN(d,k)
  seq7 <- multN(seq,k)
  di <- with(libscheme, round(probs(d7)*sum(s[-length(unique(class))])^k,0))
  pi <- probs(seq7)
  mult <- with(libscheme, s*c)
  list(data=data.frame(class = as.vector(d7), di = di, probs = pi),
       info=list(nucleotides=sum(with(libscheme, mult)), 
                 valid=with(libscheme, sum(mult[-length(mult)])),
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
#' qplot(detect(lib, size=10^8), weight=di, geom="histogram", data=lib$data)
detect <- function(lib = NNK(7), size = 10^8) {
  with(lib$data, 1 - exp(-size*probs/di))
}


getNeighborOne <- function(x, blosum=1) {
  ## For CRAN check
  BLOSUM80 <- AA1 <- Blosum <- AA2 <- NULL
    
  data(BLOSUM80, envir=environment())
  
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
  ## For CRAN check
  BLOSUM80 <- AA1 <- Blosum <- AA2 <- NULL
    
  data(BLOSUM80, envir=environment())
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
#' Use this function for only a few peptide sequences. Any larger number of peptide sequences will take too much main memory.
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
  data(BLOSUM80, envir=environment())
  if (length(x) == 1) return(getNofNeighborsOne(x, blosum, method, libscheme))

  return(llply(x, getNofNeighborsOne, blosum, method, libscheme))
}




#' Compute the number of codons for a vector of peptide sequences
#' 
#' use this function for only a few peptide sequences. Any larger number of peptide sequences should be dealt with in the framework of the library scheme and the detect function.
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
#' use this function for only a few peptide sequences. Any larger number of peptide sequences should be dealt with in the framework of the library scheme and the detect function.
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
NULL