\name{fdisp_w}
\alias{fdisp_w}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fdisp_w(d, a, tol = 1e-07, Weigthedby = c("abundance", "biomasCarabids", "biomasBees", "biomassValue"), biomassValue = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
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
function (d, a, tol = 1e-07, Weigthedby = c("abundance", "biomasCarabids", 
    "biomasBees", "biomassValue"), biomassValue = NA) 
{
    if (!inherits(d, "dist")) 
        stop("'d' must be a 'dist' object.")
    n <- attr(d, "Size")
    if (is.null(attr(d, "Labels"))) 
        stop("'d' must have labels.", "\n")
    else sn.d <- attr(d, "Labels")
    if (missing(a)) {
        ab.names <- list("Community1", sn.d)
        a <- matrix(1, 1, n, dimnames = ab.names)
    }
    com <- nrow(a)
    if (!is.matrix(a)) 
        stop("'a' must be a matrix.")
    if (ncol(a) != n) 
        stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
    if (is.null(colnames(a))) 
        stop("'a' must have column names", "\n")
    else sn.a <- colnames(a)
    if (any(sn.d != sn.a)) 
        stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
            "\n")
    a[which(is.na(a))] <- 0
    abun.sum <- apply(a, 1, sum)
    if (any(abun.sum == 0)) 
        stop("At least one community has zero-sum abundances (no species).", 
            "\n")
    abun.sum2 <- apply(a, 2, sum)
    if (any(abun.sum2 == 0)) 
        stop("At least one species does not occur in any community (zero total abundance across all communities).", 
            "\n")
    if (any(is.na(d))) 
        stop("NA's in the distance matrix.", "\n")
    A <- matrix(0, ncol = n, nrow = n)
    A[row(A) > col(A)] <- -0.5 * d^2
    A <- A + t(A)
    G <- bicenter.wt(A)
    e <- eigen(G, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    w0 <- eig[n]/eig[1]
    if (w0 > -tol) 
        r <- sum(eig > (eig[1] * tol))
    else r <- length(eig)
    vectors <- vectors[, 1:r, drop = FALSE] \%*\% diag(sqrt(abs(eig <- eig[1:r])), 
        r)
    dimnames(vectors) <- list(colnames(a), NULL)
    pos <- eig > 0
    avg.dist.cent <- rep(NA, nrow(a))
    names(avg.dist.cent) <- row.names(a)
    for (i in 1:com) {
        pres <- which(a[i, ] > 0)
        nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
        if (nb.sp >= 2) {
            w <- a[i, pres]
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
                if (is.vector(biomassValue2) == FALSE) {
                  w <- w * biomassValue2[i, pres]
                }
                else {
                  w <- w * biomassValue2[pres]
                }
            }
            centroid <- apply(vec, 2, weighted.mean, w = w)
            dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
            dist.pos <- rowSums(dist.pos^2)
            if (any(!pos)) {
                dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
                dist.neg <- rowSums(dist.neg^2)
            }
            else dist.neg <- 0
            zij <- sqrt(abs(dist.pos - dist.neg))
            avg.dist.cent[i] <- weighted.mean(zij, w)
        }
        else avg.dist.cent[i] <- 0
    }
    return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
