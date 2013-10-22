\name{FD_Clark}
\alias{FD_Clark}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FD_Clark(S, A, w = NA, Distance.method = "gower", ord = c("podani", "metric"), Cluster.method = c(ward = "ward", single = "single", complete = "complete", UPGMA = "average", UPGMC = "centroid", WPGMC = "median", WPGMA = "mcquitty"), stand.x = TRUE, stand.FD = FALSE, Weigthedby = c("abundance", "biomasCarabids", "biomasBees", "biomassValue"), biomassValue = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{Distance.method}{
%%     ~~Describe \code{Distance.method} here~~
}
  \item{ord}{
%%     ~~Describe \code{ord} here~~
}
  \item{Cluster.method}{
%%     ~~Describe \code{Cluster.method} here~~
}
  \item{stand.x}{
%%     ~~Describe \code{stand.x} here~~
}
  \item{stand.FD}{
%%     ~~Describe \code{stand.FD} here~~
}
  \item{Weigthedby}{
%%     ~~Describe \code{Weigthedby} here~~
}
  \item{biomassValue}{
%%     ~~Describe \code{biomassValue} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (S, A, w = NA, Distance.method = "gower", ord = c("podani", 
    "metric"), Cluster.method = c(ward = "ward", single = "single", 
    complete = "complete", UPGMA = "average", UPGMC = "centroid", 
    WPGMC = "median", WPGMA = "mcquitty"), stand.x = TRUE, stand.FD = FALSE, 
    Weigthedby = c("abundance", "biomasCarabids", "biomasBees", 
        "biomassValue"), biomassValue = NA) 
{
    require(FD)
    require(cluster)
    require(vegan)
    Out <- data.frame(comm = rep(NA, nrow(A)), n_sp = rep(NA, 
        nrow(A)), n_tr = rep(NA, nrow(A)), qual.FD = rep(NA, 
        nrow(A)), FDabund = rep(NA, nrow(A)), FDjointabund = rep(NA, 
        nrow(A)))
    Out$comm <- rownames(A)
    Out$n_tr <- ncol(S)
    if (is.na(w)[1]) {
        w <- rep(1, ncol(S))
    }
    if (Distance.method == "gower") {
        D <- gowdis(S, w = w, ord = ord)
    }
    else {
        if (stand.x == TRUE) {
            S2 <- scale(S, center = TRUE, scale = TRUE)
            D <- dist(S2, method = Distance.method)
        }
        else {
            D <- dist(S, method = Distance.method)
        }
    }
    tree <- hclust(D, method = Cluster.method)
    plot(tree)
    xtree <- Xtree(tree)
    c_distance <- cophenetic(tree)
    Out[, 4] <- rep(cor(D, c_distance), nrow(Out))
    if (Weigthedby != "abundance") {
        if (Weigthedby == "biomasCarabids") {
            biomassValue2 <- Jelaska(biomassValue)
        }
        if (Weigthedby == "biomasBees") {
            biomassValue2 <- Cane(biomassValue)
        }
        else {
            biomassValue2 <- biomassValue
        }
        AA <- A * biomassValue2
    }
    else AA <- A
    AFw <- AA
    for (k in 1:nrow(AA)) {
        AFw[k, ] <- AA[k, ]/max(AA[k, ])
    }
    for (i in 1:nrow(A)) {
        Out[i, 2] <- length(A[i, -which(A[i, ] == 0)])
        if (Distance.method != "gower" & stand.x == TRUE) {
            S2 <- scale(S, center = TRUE, scale = TRUE)
        }
        else {
            S2 <- S
        }
        num <- sapply(S2, is.numeric)
        Snum <- S2[, num]
        Sabund <- Snum
        if (stand.x == TRUE) {
            Sabund <- scale(Sabund, center = TRUE, scale = TRUE)
        }
        if (is.data.frame(Snum) == TRUE) {
            for (k in 1:ncol(Snum)) {
                Sabund[, k] <- Snum[, k] * as.numeric(AFw[i, 
                  ])
            }
        }
        else {
            for (k in 1:length(Snum)) {
                Sabund <- Snum[k] * as.numeric(AFw[i, ])
            }
        }
        Sabund <- cbind(Sabund, S2[, !num])
        Sabund <- Sabund[which(A[i, ] > 0), ]
        wabund <- rep(1, ncol(Sabund))
        if (Distance.method == "gower") {
            is.bin <- function(x) all(x[!is.na(x)] \%in\% c(0, 
                1))
            bin.var <- rep(NA, dim(Sabund)[2])
            names(bin.var) <- colnames(Sabund)
            for (j in 1:dim(Sabund)[2]) bin.var[j] <- is.bin(Sabund[, 
                j])
            if (any(bin.var == TRUE)) {
                Sabund[, which(bin.var == TRUE)] <- sapply(Sabund[, 
                  which(bin.var == TRUE)], as.numeric)
            }
            else {
                Dabund <- gowdis(Sabund, w = wabund, ord = ord)
            }
        }
        else {
            Dabund <- dist(Sabund, method = Distance.method)
        }
        if (attr(Dabund, "Size") <= 2) {
            print("FDabund not calculated for <= 2 species. NA inserted")
            Out[i, 5] <- NA
        }
        else {
            treeabund <- hclust(Dabund, method = Cluster.method)
            xtreeabund <- Xtree(treeabund)
            Acomm <- A[i, which(A[i, ] > 0)]
            species_in_Cabund <- ifelse(Acomm > 0, 1, 0)
            i.primeCabund <- ifelse(colSums(xtreeabund$H1[which(species_in_Cabund > 
                0), ]) > 0, 1, 0)
            Out[i, 5] <- sum(i.primeCabund * xtreeabund$h2.prime)
        }
        Dw <- D
        D2 <- as.matrix(D)
        Dw2 <- as.matrix(Dw)
        for (a in 1:ncol(AFw)) {
            for (j in 1:ncol(AFw)) {
                Dw2[a, j] <- D2[a, j] * (AFw[i, a] * AFw[i, j])
            }
        }
        sp_present <- which(A[i, ] > 0)
        Dw2 <- Dw2[sp_present, sp_present]
        Dw <- as.dist(Dw2)
        if (attr(Dw, "Size") <= 2) {
            print("FDjointabund not calculated for <= 2 species. NA inserted")
            Out[i, 6] <- NA
        }
        else {
            tree_jointabund <- hclust(Dw, method = Cluster.method)
            xtree_jointabund <- Xtree(tree_jointabund)
            species_in_Cjointabund <- ifelse(Acomm > 0, 1, 0)
            i.primeC_jointabund <- ifelse(colSums(xtree_jointabund$H1[which(species_in_Cjointabund > 
                0), ]) > 0, 1, 0)
            Out[i, 6] <- sum(i.primeC_jointabund * xtree_jointabund$h2.prime)
        }
    }
    Out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
