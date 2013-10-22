\name{FD_dendro}
\alias{FD_dendro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FD_dendro(S, A, w = NA, Distance.method = "gower", ord = c("podani", "metric"), Cluster.method = c(ward = "ward", single = "single", complete = "complete", UPGMA = "average", UPGMC = "centroid", WPGMC = "median", WPGMA = "mcquitty"), stand.x = TRUE, stand.FD = FALSE, Weigthedby = c("abundance", "biomasCarabids", "biomasBees", "biomassValue"), biomassValue = NA)
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
        nrow(A)), n_tr = rep(NA, nrow(A)), FDpg = rep(NA, nrow(A)), 
        FDw = rep(NA, nrow(A)), FDwcomm = rep(NA, nrow(A)), qual.FD = rep(NA, 
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
    Out[, 7] <- rep(cor(D, c_distance), nrow(Out))
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
    AFcomm <- AA
    for (k in 1:nrow(AA)) {
        AFcomm[k, ] <- AA[k, ]/max(AA)
    }
    for (i in 1:nrow(A)) {
        species_in_C <- ifelse(A[i, ] > 0, 1, 0)
        select_xtree <- xtree$H1[which(species_in_C > 0), ]
        if (is.array(select_xtree) == TRUE) {
            i.primeC <- ifelse(colSums(select_xtree) > 0, 1, 
                0)
        }
        else {
            i.primeC <- select_xtree
        }
        Out[i, 4] <- sum(i.primeC * xtree$h2.prime)
        Out[i, 2] <- length(A[i, -which(A[i, ] == 0)])
        xtree.weigths <- xtree$H1
        for (k in 1:nrow(S)) {
            xtree.weigths[k, ] <- ifelse(xtree$H1[k, ] > 0, AFw[i, 
                k], 0)
        }
        i.primeW <- c(1:ncol(xtree.weigths))
        for (k in 1:ncol(xtree.weigths)) {
            if (sum(xtree.weigths[which(xtree.weigths[, k] > 
                0), k]) != 0) {
                i.primeW[k] <- prod(xtree.weigths[which(xtree.weigths[, 
                  k] > 0), k], na.rm = TRUE)
            }
            else {
                i.primeW[k] <- 0
            }
        }
        Out[i, 5] <- sum(i.primeW * xtree$h2.prime)
        xtree.weigths <- xtree$H1
        for (k in 1:nrow(S)) {
            xtree.weigths[k, ] <- ifelse(xtree$H1[k, ] > 0, AFcomm[i, 
                k], 0)
        }
        i.primeW <- c(1:ncol(xtree.weigths))
        for (k in 1:ncol(xtree.weigths)) {
            if (sum(xtree.weigths[which(xtree.weigths[, k] > 
                0), k]) != 0) {
                i.primeW[k] <- prod(xtree.weigths[which(xtree.weigths[, 
                  k] > 0), k], na.rm = TRUE)
            }
            else {
                i.primeW[k] <- 0
            }
        }
        Out[i, 6] <- sum(i.primeW * xtree$h2.prime)
    }
    if (stand.FD == TRUE) {
        Out[, 4] <- Out[, 4]/max(Out[, 4])
    }
    Out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
