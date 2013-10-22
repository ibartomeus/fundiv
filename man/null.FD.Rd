\name{null.FD}
\alias{null.FD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
null.FD(S, A, it, w = NA)
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
  \item{w}{
%%     ~~Describe \code{w} here~~
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
function (S, A, it, w = NA) 
{
    require(GLDEX)
    rich_lev <- rowSums(ifelse(A > 0, 1, 0))
    rich_loop <- unique(rich_lev)
    A2 <- A
    for (r in rich_loop) {
        null_com <- matrix(ncol = ncol(A), nrow = it, data = 0)
        colnames(null_com) <- colnames(A)
        rownames(null_com) <- rep(paste("null", r, sep = "_"), 
            it)
        for (i in 1:it) {
            null_com[i, c(sample(c(1:ncol(A)), r, replace = FALSE))] <- 1
            for (j in 1:ncol(A)) {
                null_com[i, j] <- ifelse(null_com[i, j] == 1, 
                  sample(fun.zero.omit(A[, j]), 1), null_com[i, 
                    j])
            }
        }
        A2 <- rbind(A2, null_com)
    }
    fun <- FDi(S = S, A = A2, Distance.method = "gower", ord = "podani", 
        w = w, stand.FRic = TRUE, stand.FD = TRUE, corr = "cailliez")
    pFD <- c()
    null_meanFD <- c()
    null_sdFD <- c()
    for (k in 1:nrow(A)) {
        null_meanFD[k] <- mean(fun[tail(which(fun$n_sp == fun$n_sp[k]), 
            it), 4])
        null_sdFD[k] <- sd(fun[tail(which(fun$n_sp == fun$n_sp[k]), 
            it), 4])
        pFD[k] <- pnorm(fun[k, 4], mean = mean(fun[tail(which(fun$n_sp == 
            fun$n_sp[k]), it), 4]), sd = sd(fun[tail(which(fun$n_sp == 
            fun$n_sp[k]), it), 4]))
    }
    pFrich <- c()
    null_meanFrich <- c()
    null_sdFrich <- c()
    for (k in 1:nrow(A)) {
        null_meanFrich[k] <- mean(fun[tail(which(fun$n_sp == 
            fun$n_sp[k]), it), 9])
        null_sdFrich[k] <- sd(fun[tail(which(fun$n_sp == fun$n_sp[k]), 
            it), 9])
        pFrich[k] <- pnorm(fun[k, 8], mean = mean(fun[tail(which(fun$n_sp == 
            fun$n_sp[k]), it), 9]), sd = sd(fun[tail(which(fun$n_sp == 
            fun$n_sp[k]), it), 9]))
    }
    out <- data.frame(comm = fun$comm[1:nrow(A)], Rich = fun$n_sp[1:nrow(A)], 
        FD = fun$FD_PG[1:nrow(A)], null_meanFD = null_meanFD, 
        null_sdFD = null_sdFD, pFD = round(pFD, 4), Frich = fun$Frich[1:nrow(A)], 
        null_meanFrich = null_meanFrich, null_sdFrich = null_sdFrich, 
        pFD = round(pFrich, 4))
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
