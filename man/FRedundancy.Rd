\name{FRedundancy}
\alias{FRedundancy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FRedundancy(S, A, it = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{it}{
%%     ~~Describe \code{it} here~~
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
function (S, A, it = 10) 
{
    rob <- data.frame(site = rownames(A), rob_FDpg = rep(NA, 
        nrow(A)), rob_FRich = rep(NA, nrow(A)))
    for (s in 1:nrow(A)) {
        Ano0 <- A[s, -which(A[s, ] == 0)]
        l_funFRich <- matrix(ncol = length(Ano0), nrow = it)
        l_funFD_PG <- matrix(ncol = length(Ano0), nrow = it)
        for (k in 1:it) {
            Ano0 <- A[s, -which(A[s, ] == 0)]
            Ar <- matrix(nrow = length(Ano0), ncol = nrow(S), 
                dimnames = list(paste("d", seq(0, length(Ano0) - 
                  1, 1), sep = ""), rownames(S)))
            Ar[1, ] <- as.numeric(A[s, ])
            for (i in seq(2, length(Ano0) - 1, 1)) {
                Ar[i, ] <- Ar[i - 1, ]
                Ar[i, which(colnames(Ar) == names(sample(Ano0, 
                  1)))] <- 0
                Ano0 <- Ar[i, -which(Ar[i, ] == 0)]
            }
            Ar[i + 1, ] <- rep(1, ncol(Ar))
            fun <- FDi(S = S, A = Ar, Distance.method = "gower", 
                ord = "podani", w = w, stand.FRic = TRUE, corr = "cailliez")
            fun[i, 8] <- 0
            l_funFRich[k, ] <- fun[, 8]
            l_funFD_PG[k, ] <- fun[, 4]
        }
        MeanFunFRich <- colMeans(l_funFRich)[-(i + 1)]
        MeanFunFD_PG <- colMeans(l_funFD_PG)[-(i + 1)]
        functFRich <- c()
        functFD_PG <- c()
        for (a in seq(from = nrow(fun) - 1, to = 2, by = -1)) {
            functFRich[which(a == seq(from = nrow(fun) - 1, to = 2, 
                by = -1))] <- abs(MeanFunFRich[a] - MeanFunFRich[a - 
                1])
            functFD_PG[which(a == seq(from = nrow(fun) - 1, to = 2, 
                by = -1))] <- abs(MeanFunFD_PG[a] - MeanFunFD_PG[a - 
                1])
        }
        functFRich[i] <- 0
        functFD_PG[i] <- 0
        ex.FD <- data.frame(no = c(1:(nrow(fun) - 1)), ext.lower = rep(1, 
            nrow(fun) - 1), funct = functFD_PG)
        ex.FD <- as.matrix(ex.FD)
        class(ex.FD) <- "bipartite"
        attr(ex.FD, "exterminated") <- "lower"
        ex.FRich <- data.frame(no = c(1:(nrow(fun) - 1)), ext.lower = rep(1, 
            nrow(fun) - 1), funct = functFRich)
        ex.FRich <- as.matrix(ex.FRich)
        class(ex.FRich) <- "bipartite"
        attr(ex.FRich, "exterminated") <- "lower"
        rob[s, 2] <- robustness(ex.FD)
        rob[s, 3] <- robustness(ex.FRich)
    }
    rob
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
