library(peptider)
data(BLOSUM80)
for (lib in c("NNKC", "NNNC", "TrimerC", "NNBC")) {

    lib <- "NNBC"
k <- 10

sch <- scheme(lib)

schl <- unlist(strsplit(as.character(sch$aacid[-length(sch$aacid)]), ""))
schd <- data.frame(AA=schl, C0=codons(schl, lib=lib))
schd$C1 <- getNofNeighbors(schd$AA, method="codon", lib=lib) - schd$C0
schd$C1T0 <- with(schd, C1/C0)

ctabs <- function(values, labels) {
    labels <- as.character(labels)
    
    tv <- xtabs(~values)
    require(plyr)
    ct <- ldply(names(tv), function(x) {
        data.frame(AA=paste(labels[which(values == x)], collapse=""),
                   c=x)
    })
    ct
}

user <- with(schd, ctabs(paste(C0, C1, sep=":"), AA))
user <- data.frame(user, ldply(strsplit(as.character(user$c), ":"), 
                               function(x) as.numeric(unlist(x))))
names(user)[3:4] <- c("c0", "c1")
user$cr <- with(user, c1/c0)
user$c <- user$c0
user$s <- nchar(as.character(user$AA))
user$one <- 1

i <- 1
dx <- user
dx$o <- dx$one
dx$L <- user$AA
while (i < k) {
    C0 <- as.vector(outer(dx$c0, user$c, FUN = "*"))
    S <- as.vector(outer(dx$s, user$s, FUN = "*"))
    CR <- as.vector(outer(dx$cr, user$cr, FUN = "+"))
    O <- as.vector(outer(dx$o, user$one, FUN = "*"))
    #    L <- as.vector(outer(L, user$AA, FUN = function(x,y) {        
    #        x <- as.character(x)
    #        y <- as.character(y)
    #       paste(pmin(x,y), pmax(x,y), sep=",") 
    #    }))    
    L <- as.vector(outer(dx$L, user$AA, "paste", sep = ","))
    lx <- strsplit(as.character(L), split=",")
    L <- laply(llply(lx, sort), paste, collapse=",")
    
    x <- data.frame(L, C0,CR,S,O)
    dx <- ddply(x, .(L), summarise, 
                c0=C0[1],
                cr=CR[1],
                s=S[1],
                o=sum(O))
    
    #    save(L,C0,CR,S,O, file=sprintf("%s-%d.RData", lib, i))
    dx$N1 <- with(dx, c0*(cr+1))
    dx$lib <- lib
    dx$k <- i+1
    cat(sprintf("stage %d done\n",i))
    write.table(dx, file="neighbors2.csv", sep=",", col.names=!file.exists("neighbors2.csv"), append=TRUE, row.names=FALSE)
    
    i <- i + 1
}
}


############

nb <- read.csv("neighbors2.csv")
save(nb, "nb.RData")

