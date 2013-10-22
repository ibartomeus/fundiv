\name{Xtree}
\alias{Xtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Xtree(h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
%%     ~~Describe \code{h} here~~
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
function (h) 
{
    species.names <- h$labels
    H1 <- matrix(0, length(h$order), 2 * length(h$order) - 2)
    l <- vector("numeric", 2 * length(h$order) - 2)
    for (i in 1:(length(h$order) - 1)) {
        if (h$merge[i, 1] < 0) {
            l[2 * i - 1] <- h$height[order(h$height)[i]]
            H1[-h$merge[i, 1], 2 * i - 1] <- 1
        }
        else {
            l[2 * i - 1] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 
                1]]]
            H1[, 2 * i - 1] <- H1[, 2 * h$merge[i, 1] - 1] + 
                H1[, 2 * h$merge[i, 1]]
        }
        if (h$merge[i, 2] < 0) {
            l[2 * i] <- h$height[order(h$height)[i]]
            H1[-h$merge[i, 2], 2 * i] <- 1
        }
        else {
            l[2 * i] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 
                2]]]
            H1[, 2 * i] <- H1[, 2 * h$merge[i, 2] - 1] + H1[, 
                2 * h$merge[i, 2]]
        }
    }
    dimnames(H1) <- list(species.names, NULL)
    list(h2.prime = l, H1 = H1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
