library(peptider)
data(BLOSUM80)

lib <- "NNKC"
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

i <- 1
C0 <- user$c
CR <- user$cr
L <- user$AA
S <- user$s
while (i < k) {
    C0 <- outer(C0, user$c, FUN = "*")
    S <- outer(S, user$s, FUN = "*")
    CR <- outer(CR, user$cr, FUN = "+")
    L <- outer(L, user$AA, FUN = "paste", sep = ",")
    
    write(C0, file=sprintf("C0-%s-%d.Rdata", lib, i+1))
    write(S, file=sprintf("S-%s-%d.Rdata",lib, i+1))
    write(CR, file=sprintf("CR-%s-%d.Rdata",lib, i+1))
    write(L, file=sprintf("L-%s-%d.Rdata",lib, i+1))
    sprintf("stage %d done",i)
    i <- i + 1
}

for (i in 1:k) {
# read matrices back in
x <- data.frame(AA=as.vector(L), 
                c0=as.vector(C0), 
                cr=as.vector(CR), 
                s = as.vector(S))
}
## number of neighbors
x$N1 <- with(x, c0*(cr+1))
x$lib <- lib
x$k <- k

write.table(x, file="neighbors.csv", sep=",", col.names=!file.exists("neighbors.csv"), append=TRUE, row.names=FALSE)
