#' Reduce the regular encoding to an easier/faster format
#' @param class The peptide class
#' @param libscheme The scheme to use
#' @return Vector of reduced peptide encodings
encodingReduce <- function(class, libscheme) {
    sapply(class, function(x) { 
        splitted <- strsplit(x, split = "")[[1]]
        lvls <- libscheme(1)$info$scheme$class
        myTbl <- table(factor(splitted, levels = lvls[lvls != "Z"]))
        return(paste(myTbl, collapse = ","))
    })
}

#' Get the counts possible for each scheme and k
#' @param libscheme The scheme to use
#' @param k Peptide length
#' @return Character vector of possible counts for each class
getCounts <- function(libscheme, k) {
    scheme.levels <- length(libscheme(1)$info$scheme$class) - 1
    
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
#' @param Custom function representing the scheme to use
#' @return A data frame of peptide probabilities
generateCustomProbs <- function(Custom) {
    ## The Constants
    k <- 6:10
    N <- c(10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12)
    n <- as.vector(sapply(10^seq(6, 12, by = 1), `*`, seq(1.0, 9.9, by = 0.1)))
    
    lib.probs.tmp <- ldply(k, function(y) {
        df <- data.frame(Counts = getCounts(Custom, y))
        df$k <- y
        
        df
    })
    lib.probs.tmp$samp.encoding <- apply(lib.probs.tmp, 1, function(z) { paste(rep(Custom(1)$info$scheme$class, as.numeric(c(strsplit(as.character(z), split = ",")[[1]], 0))), collapse = "") })
    
    lib.probs <- ldply(k, function(y) {
        lib.data <- get("Custom")(y)$data
        lib.data$class <- gsub("\\.", "", lib.data$class)
        lib.data.subset <- subset(lib.data, class %in% lib.probs.tmp$samp.encoding)
        
        return.df <- data.frame(Counts = encodingReduce(as.character(lib.data.subset$class), Custom), samp.encoding = lib.data.subset$class, probs = lib.data.subset$probs)
        #return.df$probs <- apply(return.df, 1, function(z) {as.numeric(z[3]) / prod(getProbCorrection(z[3], x))})
        return.df$k <- y
        return.df$di <- lib.data.subset$di
        return.df$choices <- sapply(as.character(return.df$Counts), getChoices)
        
        return.df
    })
    
    return(lib.probs)
}

#' For a given scheme, generate a dataset with the library information
#' @param Custom function representing the scheme to use
#' @return A data frame of library information
generateCustomLib <- function(Custom) {
    ## The Constants
    k <- 6:10
    N <- c(10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12)
    n <- as.vector(sapply(10^seq(6, 12, by = 1), `*`, seq(1.0, 9.9, by = 0.1)))
    
    lib.stats <- ldply(k, function(k1) {
        cat("Processing for k =", k1, "\n")
        
        lib <- Custom(k1)
        n <- as.vector(sapply(10^seq(6, 12, by = 1), `*`, seq(1.0, 9.9, by = 0.1)))
        
        lib.stats <- ldply(n, .progress = "text", function(n1) {  
            # cat("Processing for n =", n1, "\n")
            
            cov = coverage(k=k1, libscheme=Custom, N=n1, lib=lib)
            eff = efficiency(k=k1, libscheme=Custom, N=n1, lib=lib)
            c(k=k1, n=n1, cov=cov, eff=eff)
        })
        
        lib.stats$div = makowski(k=k1, libscheme=Custom)
        
        lib.stats
    })
    
    return(lib.stats)
}

#' Generate peptide and library information for a given scheme
#' @param scheme_name The name of the resulting encoding scheme
#' @param scheme_def A data frame containing encoding information for the scheme
#' @return TRUE upon completion of the script and output of the CSV files
#' @export
#' @examples
#' \dontrun{
#' generateCustom()
#' generateCustom(scheme_name = "NNK_S", scheme_def = peptider:::nnk_scheme)
#' }
generateCustom <- function(scheme_name = "Custom", scheme_def = read.csv(file.choose())) {
    ## Load the scheme
    Custom <- function(k) {
        libBuild(k, scheme_def)
    }
    
    custom.probs <- generateCustomProbs(Custom)
    custom.lib <- generateCustomLib(Custom)
    
    write.csv(custom.probs, paste("prob-", scheme_name, ".csv", sep = ""), row.names = FALSE)
    write.csv(custom.lib, paste("lib-", scheme_name, ".csv", sep = ""), row.names = FALSE)
    write.csv(file, paste("scheme-", scheme_name, ".csv", sep = ""), row.names = FALSE)
    
    TRUE
}
