#' Built-in library schemes
#'
#' @name schemes
#' @title Built-in library schemes for peptider
#' @description This data set contains descriptions of amino acid classes several commonly used library schemes: NNN, NNB, NNK, trimer, and variations of each in which Cysteine is not considered a viable amino acid.
#' 
#' NNN: All four bases (\"N\" = G/A/T/C) possible at all three positions in the codon.
#' NNB: All four bases in the first two codon positions possible, the third position is restricted to G, T or C (= \"B\")
#' NNK/S: All four bases in the first two codon positions possible, the third position is restricted to G/T (= \"K\") or two C/G (= \"S\").
#' trimer: trimer describes the concept that DNA is assembled from prefabricated trimeric building blocks. This allows the generation of libraries from a predefined set of codons and thereby complete exclusion of Stop codons and other unwanted codons.
#' NNN (-C): NNN with Cysteine ignored.
#' NNB (+C): NNB with Cysteine ignored.
#' NNK/SC (+C): NNK/S with Cysteine ignored.
#' trimer (+C): Trimer with Cysteine ignored.
#' 
#' The schemes differ in the number of used codons, ranging from 64 (NNN), 48 (NNB), 32 (NNK/S) to 20 or less (trimer). Coding schemes that allow varying ratios of codons/amino acid, result in libraries biased towards amino acids which are encoded more often. Further, the number of Stop codons that can lead to premature termination of the peptide sequence influences the performance of the library.
#' @docType data
#' @usage data(schemes)
NULL

#' Get the specified library scheme definition
#' 
#' @param name name of the scheme as a character vector
#' 
#' @return list consisting of a data frame of peptide classes, amino acids, and size of the classes
#' 
#' @export
#' 
#' @examples
#' scheme("NNN")
#' scheme("NNK")
scheme <- function(name) {
    schemes <- NULL ## R Check
    data(schemes, envir=environment())
    
    scheme_def <- schemes[[paste(tolower(name), "scheme", sep = "_")]]
    if (is.null(scheme_def)) stop(paste("No library with name", name, "is included in peptider"))
    
    return(scheme_def)
}

#' Get the specified library scheme
#' 
#' @param schm either a character vector giving the name of a built-in scheme, or a data frame consisting of the scheme definition
#' @param k length of peptide sequences
#' 
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, and a list of additional information relating to the library scheme
#' 
#' @export
#' 
#' @examples
#' libscheme("NNN")
#' libscheme("NNK", 2)
#' 
#' # Build a custom trimer library
#' custom <- data.frame(class = c("A", "Z"), aacids = c("SLRAGPTVIDEFHKNQYMW", "*"), c = c(1, 0))
#' libscheme(custom)
libscheme <- function(schm, k = 1) {
    if (is.character(schm)) return(libBuild(scheme(schm), k = k))
    else if (is.data.frame(schm)) return(libBuild(k, schm))
    else stop("scheme must be either a character or a data frame")
}

#' Diversity index according to Makowski
#'
#' Diversity according to Makowski is defined as ... need a reference to Makowski & Soares 2003 here. 
#' @param k length of peptide sequences
#' @param libscheme Name (character vector) or definition (data frame) of scheme
#' @return diversity index between 0 and 1
#' @export
#' @examples
#' makowski(2, "NNN")
#' makowski(3, "NNK")
#' makowski(3, "Trimer")
makowski <- function(k, libscheme) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    scheme_def <- libscheme(libschm, k)
    
    dframe <- scheme_def$data
    info <- scheme_def$info$scheme
    numAA <- sum(info$s[-nrow(info)]) 
    
    with(dframe, 1/(numAA^k*sum(probs^2/di)))
}

#' Coverage as expected number of peptides given all possible peptides
#'
#' Coverage of library of size N given random sampling from the pool of all possible peptides according 
#' to probabilities determined according to the library scheme.
#' @param k length of peptide sequences
#' @param libscheme Name (character vector) or definition (data frame) of scheme
#' @param N size of the library 
#' @param lib library scheme
#' @return coverage index between 0 and 1
#' @export
#' @examples
#' coverage(2, "NNN", 10^3)
#' coverage(2, "NNK", 10^3)
#' coverage(2, "Trimer", 10^3) ## Trimer coverage is not 1 because of random sampling.
coverage <- function(k, libscheme, N, lib=NULL) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    if (is.null(lib)) lib <- libscheme(libschm, k)
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
#' @param libscheme Name (character vector) or definition (data frame) of scheme
#' @param N size of the library 
#' @param lib library, if null, libscheme will be used to create it
#' @return relative efficiency index between 0 and 1
#' @export
#' @examples
#' efficiency(3, "NNN", 10^2)
#' efficiency(3, "NNK", 10^2)
#' efficiency(3, "Trimer", 10^2) ## Trimer efficiency is not 1 because of random sampling.
efficiency <- function(k, libscheme, N, lib=NULL) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    if (is.null(lib)) lib <- libscheme(libschm, k)
    libdata <- lib$data
    
    initialloss <- (1-(lib$info$valid/lib$info$nucleotides)^k)
    libdata$expected <- libdata$probs*N*(1-initialloss)
    libdata$z <- with(libdata, di*(1-exp(-expected/di)))
    
    s_count <- sum(subset(lib$info$scheme, class != "Z")$s)
    
    with(libdata, min(s_count^k,sum(z))/N)
}

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
#' @param lib library used in experiment, defaults to NNK with peptide length 7
#' @param size size of the library, defaults to 10^8
#' @return vector of detection probabilities for peptide sequences in each class
#' @export
#' @examples
#' summary(detect())
#'
#' require(ggplot2)
#' lib <- libscheme("NNK", 7)
#' qplot(detect(lib, size=10^8), weight=di, geom="histogram", data=lib$data)
detect <- function(lib = libscheme("NNK", 7), size = 10^8) {
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
    
    stopifnot(!(is.null(libscheme) & nchar(libscheme) == 0))    
    lib <- peptider::libscheme(libscheme)$info$scheme
    
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
#' getNofNeighbors("N", method="dna", libscheme="NNK")
getNofNeighbors <- function(x, blosum = 1, method="peptide", libscheme=NULL) {
    data(BLOSUM80, envir=environment())
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    if (length(x) == 1) return(getNofNeighborsOne(x, blosum, method, libschm))
    
    return(llply(x, getNofNeighborsOne, blosum, method, libschm))
}

#' Compute the number of codons for a vector of peptide sequences
#' 
#' use this function for only a few peptide sequences. Any larger number of peptide sequences should be dealt with in the framework of the library scheme and the detect function.
#' @param x (vector) of character strings of  peptide sequences.
#' @param libscheme library scheme under which neighbors are being calculated. this is only of importance, if method="dna"
#' @param flag internal use only: Set to true if calling this from another function
#' @return vector of numbers of codons 
#' @export
#' @examples
#' codons("APE", libscheme="NNK")
#' codons("HENNING", libscheme="NNK")
codons <- function(x, libscheme=NULL, flag = FALSE) {
    libschm <- libscheme
    if (!flag) {
        libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
        if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    }
    if (length(x) == 1) return(codonsOne(x, libschm))
    
    unlist(llply(x, codonsOne, libscheme=libschm))
}

codonsOne <- function(x, libscheme) {
    stopifnot(!(is.null(libscheme) & nchar(libscheme) == 0))    
    lib <- libscheme(libscheme)$info$scheme
    
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
#' ppeptide("APE", libscheme="NNK", N=10^8)
#' ppeptide("HENNING", libscheme="NNK", N=10^8)
ppeptide <- function(x, libscheme, N) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    n <- sum(codons(x, libscheme=libschm, flag = TRUE))
    Max <- peptider::libscheme(libschm, 1)$info$valid^nchar(as.character(x[1]))
    1 - exp(-N*n/Max)
}

#' BLOSUM80 matrix
#' 
#' @name BLOSUM80
#' @title BLOSUM80 matrix
#' @description The BLOSUM80 matrix, which stands for Blocks Substitution Matrix, defines log-odds scores for the ratio of the chance of two amino acids appearing in a sequence over the chance that the two amino acids appear in any sequence.  Larger scores indicate a higher probability of substitutions.  This matrix is used in order to compute sequences which are in the neighborhood of other sequences.
#' @docType data
#' @usage data(BLOSUM80)
NULL