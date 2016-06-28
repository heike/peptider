library(plyr)
library(peptider)

leading.zeros <- function(vec) {
    count <- 0
    for (i in 1:length(vec)) {
        val <- vec[i]
        if (val == 0) count <- count + 1 else break 
    }
    
    return(count)
}

scheme.list <- c("NNN" = "NNN", "NNB" = "NNB", "NNK/S" = "NNK", "trimer" = "trimer", "NNN (-C)" = "NNNC", "NNB (-C)" = "NNBC", "NNK/S (-C)" = "NNKC", "trimer (-C)" = "trimerC")
k <- 1:20
N <- as.vector(sapply(10^seq(1, 25, by = 1), `*`, seq(1.0, 9.9, by = 0.1)))

test_routine <- function(lib, k, scheme) {   
    total_N <- as.vector(sapply(10^seq(1, 25, by = 1), `*`, seq(1.0, 9.9, by = 0.1)))
    
    s_count <- sum(lib$info$scheme$s[lib$info$scheme$class != "Z"])
    
    mydivs <- sapply(10^(1:25), function(y){peptider::diversity(k = k, libscheme = scheme, N = y, lib = lib)})
    covs <- mydivs / s_count^k
    effs <- pmin(mydivs / 10^(1:25), 1)
    
    new_N <- as.vector(sapply(10^seq((which.max(effs)), 25, by = 1), `*`, seq(1.0, 9.9, by = 0.1)))
    mydivs <- sapply(new_N, function(y){peptider::diversity(k = k, libscheme = scheme, N = y, lib = lib)})
    covs <- mydivs / s_count^k
    effs <- pmin(mydivs / new_N, 1)
    effs[1:(length(effs) - which.max(rev(effs)))] <- max(effs)
    
    if (length(mydivs) < length(total_N)) covs <- c(rep(0, length(total_N) - length(mydivs)), covs)
    if (length(mydivs) < length(total_N)) effs <- c(rep(max(effs), length(total_N) - length(mydivs)), effs)
    
    data.frame(scheme = scheme, k = k, N = total_N, coverage = covs, efficiency = effs)
}

generateLibData <- function(scheme, k, N) {
    dframe <- expand.grid(scheme = scheme, k = k)
    dframe$scheme <- as.character(dframe$scheme)
    
    libframe <- ddply(dframe, .(scheme, k), .progress = "text", function(x) {
        mylib <- peptider::libscheme(x$scheme, x$k)
        
        test_routine(mylib, x$k, x$scheme)
    })
    
    return(libframe)
}

generatePepData <- function(scheme, k) {
    dframe <- expand.grid(scheme = scheme, k = k)
    dframe$scheme <- as.character(dframe$scheme)
    
    pepframe <- ddply(dframe, .(scheme, k), function(x) {
        mylib <- peptider:::libscheme(x$scheme, x$k)
        mylibdata <- mylib$data
        
        mylibdata
    })
    
    pepframe
}

generateNeiData <- function(scheme, k) {
    dframe <- expand.grid(scheme = scheme, k = k)
    dframe$scheme <- as.character(dframe$scheme)
    
    nei.dat <- ddply(dframe, .(scheme, k), function(x) {
            
        my.schemedef <- scheme(x$scheme)
        
        worst.aa <- strsplit(as.character(my.schemedef$aacid)[which.min(my.schemedef$c[my.schemedef$c > 0])], "")[[1]]
        best.aa <- strsplit(as.character(my.schemedef$aacid)[which.max(my.schemedef$c)], "")[[1]]
        
        peps.worst <- sapply(worst.aa, function(aa){paste(rep(aa, x$k), collapse = "")})
        peps.best <- sapply(best.aa, function(aa){paste(rep(aa, x$k), collapse = "")})
        
        worst.tmp <- sapply(peps.worst, function(pep){
            ppeptide(getNeighbors(pep), x$scheme, N)
        })
        worst <- apply(worst.tmp, 1, min)
        worst.which <- which.min(worst.tmp[1,])
        worst2 <- ppeptide(unique(unlist(getNeighbors(getNeighbors(peps.worst[worst.which])))), x$scheme, N)
        
        best.tmp <- sapply(peps.best, function(pep){
            ppeptide(getNeighbors(pep), x$scheme, N)
        })
        best <- apply(best.tmp, 1, max)
        best.which <- which.max(best.tmp[1,])
        best2 <- ppeptide(unique(unlist(getNeighbors(getNeighbors(peps.best[best.which])))), x$scheme, N)
        
        data.frame(scheme = x$scheme, k = x$k, N = rep(N, times = 2), degree = rep(1:2, each = length(N)), worst = c(worst, worst2), best = c(best, best2))
    })
    
    nei.dat
}

library_data <- generateLibData(as.character(scheme.list), k, N)
peptide_data <- generatePepData(as.character(scheme.list), k)
neighborhood_data <- generateNeiData(as.character(scheme.list), k)

write.csv(library_data, "data/library_data.csv", row.names = FALSE)
write.csv(peptide_data, "data/peptide_data.csv", row.names = FALSE)
write.csv(neighborhood_data, "data/neighborhood_data.csv", row.names = FALSE)
