#' Reduce the regular encoding to an easier/faster format
#' @param class The peptide class
#' @param libscheme The scheme to use
#' @return Vector of reduced peptide encodings
encodingReduce <- function(class, libscheme) {
    sapply(class, function(x) { 
        splitted <- strsplit(x, split = "")[[1]]
        lvls <- libscheme$info$scheme$class
        myTbl <- table(factor(splitted, levels = lvls[lvls != "Z"]))
        return(paste(myTbl, collapse = ","))
    })
}

#' Get the counts possible for each scheme and k
#' @param libscheme The scheme to use
#' @param k Peptide length
#' @return Character vector of possible counts for each class
getCounts <- function(libscheme, k) {
    scheme.levels <- length(libscheme$info$scheme$class) - 1
    
    str.func <- paste("expand.grid(", paste(rep(paste("0:", k, sep = ""), times = scheme.levels), collapse = ", "), ")")
    
    grid.df <- eval(parse(text = str.func))
    grid.sub <- subset(grid.df, apply(grid.df, 1, sum) == k)
    grid.str <- apply(grid.sub, 1, paste, collapse = ",")
    
    return(as.character(grid.str))
}

#' Get the number of peptides that reduce to a particular reduced encoding
#' @param str The reduced encoding string
#' @return An integer of the possible number of peptides reducing to this encoding
getChoices <- function(str) {
    test <- as.numeric(unlist(strsplit(str, split = ",")))
    
    left <- sum(test)
    total <- 1
    for (i in 1:length(test)) {
        val <- choose(left, test[i])
        left <- left - test[i]
        total <- total * val
    }
    
    return(total)
}

#' For a given scheme, generate a dataset with the peptide probabilities
#' @param scheme_def definition of the custom scheme
#' @param k peptide lengths to include
#' @param n exponents of the library size to include
#' @import plyr
#' @return A data frame of peptide probabilities
generateCustomProbs <- function(scheme_def, k = 6:10, n = 6:14) {
    ## Library sizes
    n <- as.vector(sapply(10^n, `*`, seq(1.0, 9.9, by = 0.1)))
    
    ## Generate scheme
    lib <- libscheme(scheme_def)
    
    cat("Getting possible peptide encodings...\n")
    lib.probs.tmp <- ldply(k, function(y) {
        df <- data.frame(Counts = getCounts(lib, y))
        df$k <- y
        
        df
    })
    
    cat("Getting a sample peptide encoding...\n")
    lib.probs.tmp$samp.encoding <- apply(lib.probs.tmp, 1, function(z) { paste(rep(lib$info$scheme$class, as.numeric(c(strsplit(as.character(z), split = ",")[[1]], 0))), collapse = "") })
    
    cat("Processing probabilities...\n")
    lib.probs <- ldply(k, .progress = "text", .fun = function(y) {
        lib.data <- libscheme(scheme_def, y)$data
        lib.data$class <- gsub("\\.", "", lib.data$class)
        lib.data.subset <- subset(lib.data, class %in% lib.probs.tmp$samp.encoding)
        
        return.df <- data.frame(Counts = encodingReduce(as.character(lib.data.subset$class), lib), samp.encoding = lib.data.subset$class, probs = lib.data.subset$probs)
        #return.df$probs <- apply(return.df, 1, function(z) {as.numeric(z[3]) / prod(getProbCorrection(z[3], x))})
        return.df$k <- y
        return.df$di <- lib.data.subset$di
        return.df$choices <- sapply(as.character(return.df$Counts), getChoices)
        
        return.df
    })
    
    return(lib.probs)
}

generateCustomProbs_new <- function(scheme_def, k = 6:10) {
    cat("Getting possible peptide encodings...\n")
    lib.probs.tmp <- ldply(k, function(y) {
        ## Generate scheme
        df <- libscheme_new(scheme_def, y)$data
        df$k <- y
        
        df
    })
    lib.probs.tmp$scheme <- "Custom"

    return(lib.probs.tmp)
}

#' For a given scheme, generate a dataset with the library information
#' @param scheme_def definition of the custom scheme
#' @param k peptide lengths to include
#' @param n exponents of the library size to include
#' @import plyr
#' @return A data frame of library information
generateCustomLib <- function(scheme_def, k = 6:10, n = 6:14) {
    ## Library sizes
    n <- as.vector(sapply(10^n, `*`, seq(1.0, 9.9, by = 0.1)))
    
    cat("Generating library properties...\n")
    lib.stats <- ldply(k, function(k1) {
        cat("Processing for k =", k1, "\n")
        
        lib <- libscheme(scheme_def, k = k1)
        
        lib.stats <- ldply(n, .progress = "text", function(n1) {  
            # cat("Processing for n =", n1, "\n")
            
            cov = coverage(k=k1, libscheme=scheme_def, N=n1, lib=lib)
            eff = efficiency(k=k1, libscheme=scheme_def, N=n1, lib=lib)
            c(k=k1, n=n1, cov=cov, eff=eff)
        })
        
        cat("Generating library diversity...\n")
        lib.stats$div = makowski(k=k1, libscheme=scheme_def)
        
        lib.stats
    })
    
    return(lib.stats)
}

generateCustomLib_new <- function(scheme_def, k = 6:10, n = 6:14) {
    ## Library sizes
    n <- as.vector(sapply(10^n, `*`, seq(1.0, 9.9, by = 0.1)))
    
    cat("Generating library properties...\n")
    lib.stats <- ldply(k, function(k1) {
        cat("Processing for k =", k1, "\n")
        
        lib <- libscheme_new(scheme_def, k = k1)
        
        lib.stats <- ldply(n, .progress = "text", function(n1) {  
            # cat("Processing for n =", n1, "\n")
            
            cov = coverage_new(k=k1, libscheme=scheme_def, N=n1, lib=lib)
            eff = efficiency_new(k=k1, libscheme=scheme_def, N=n1, lib=lib)
            c(k=k1, N=n1, coverage=cov, efficiency=eff)
        })        
        cat("Generating library diversity...\n")
        lib.stats$diversity = makowski_new(k=k1, libscheme=scheme_def)
        lib.stats$scheme <- "Custom"
        
        lib.stats
    })
    
    return(lib.stats)
}

generateCustomNei_new <- function(scheme_def, k = 6:10, n = 6:14) {
    ## Library sizes
    n <- as.vector(sapply(10^n, `*`, seq(1.0, 9.9, by = 0.1)))
    
    nei.dat <- ldply(k, function(x) {
        
        worst.aa <- strsplit(as.character(scheme_def$aacid)[which.min(scheme_def$c)], "")[[1]]
        best.aa <- strsplit(as.character(scheme_def$aacid)[which.max(scheme_def$c)], "")[[1]]
        
        peps.worst <- sapply(worst.aa, function(aa){paste(rep(aa, x), collapse = "")})
        peps.best <- sapply(best.aa, function(aa){paste(rep(aa, x), collapse = "")})
        
        worst.tmp <- sapply(peps.worst, function(pep){
            ppeptide(getNeighbors(pep), scheme_def, n)
        })
        worst <- apply(worst.tmp, 1, min)
        worst.which <- which.min(worst.tmp[1,])
        worst2 <- ppeptide(unique(unlist(getNeighbors(getNeighbors(peps.worst[worst.which])))), scheme_def, n)
        
        best.tmp <- sapply(peps.best, function(pep){
            ppeptide(getNeighbors(pep), scheme_def, n)
        })
        best <- apply(best.tmp, 1, max)
        best.which <- which.max(best.tmp[1,])
        best2 <- ppeptide(unique(unlist(getNeighbors(getNeighbors(peps.best[best.which])))), scheme_def, n)
        
        data.frame(scheme = "Custom", k = x, N = rep(n, times = 2), degree = rep(1:2, each = length(n)), worst = c(worst, worst2), best = c(best, best2))
    })
    
    nei.dat
}

#' Generate peptide and library information for a given scheme
#' 
#' This function will generate library properties for a custom scheme.  It is primarily intended to be used on http://www.pelica.org.
#' @param scheme_name The name of the resulting encoding scheme
#' @param scheme_def A data frame containing encoding information for the scheme
#' @param k peptide lengths to include
#' @param n exponents of the library size to include
#' @return TRUE upon completion of the script and output of the CSV files
#' @export
#' @examples
#' \dontrun{
#' generateCustom()
#' generateCustom(scheme_name = "NNN", scheme_def = scheme("NNN"))
#' }
generateCustom <- function(scheme_name = "Custom", scheme_def = read.csv(file.choose()), k = 6:10, n = 6:14) {
    custom.probs <- generateCustomProbs(scheme_def, k, n)
    custom.lib <- generateCustomLib(scheme_def, k, n)
    
    write.csv(custom.probs, paste("prob-", scheme_name, ".csv", sep = ""), row.names = FALSE)
    write.csv(custom.lib, paste("lib-", scheme_name, ".csv", sep = ""), row.names = FALSE)
    write.csv(file, paste("scheme-", scheme_name, ".csv", sep = ""), row.names = FALSE)
    
    TRUE
}

generateCustom_new <- function(scheme_name = "custom", scheme_def = read.csv(file.choose()), k = 6:10, n = 6:14) {
    custom.probs <- generateCustomProbs_new(scheme_def, k)
    custom.lib <- generateCustomLib_new(scheme_def, k, n)
    custom.nei <- generateCustomNei_new(scheme_def, k, n)
    
    write.csv(custom.probs, paste("peptide_", scheme_name, ".csv", sep = ""), row.names = FALSE)
    write.csv(custom.lib, paste("library_", scheme_name, ".csv", sep = ""), row.names = FALSE)
    write.csv(custom.nei, paste("neighborhood_", scheme_name, ".csv", sep = ""), row.names = FALSE)
    write.csv(scheme_def, paste("scheme_", scheme_name, ".csv", sep = ""), row.names = FALSE)
    
    TRUE
}
