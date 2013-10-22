\name{functcomp_bm}
\alias{functcomp_bm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
functcomp_bm(x, a, CWM.type = c("dom", "all"), bin.num = NULL, Weigthedby = c("abundance", "biomasCarabids", "biomasBees", "biomassValue"), biomassValue = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{CWM.type}{
%%     ~~Describe \code{CWM.type} here~~
}
  \item{bin.num}{
%%     ~~Describe \code{bin.num} here~~
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
function (x, a, CWM.type = c("dom", "all"), bin.num = NULL, Weigthedby = c("abundance", 
    "biomasCarabids", "biomasBees", "biomassValue"), biomassValue = NA) 
{
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("'x' must be a matrix or a data frame.", "\n")
    else x <- data.frame(x)
    if (!is.matrix(a)) 
        stop("'a' must be a matrix.", "\n")
    if (is.null(row.names(x))) 
        stop("'x' must have row names.", "\n")
    else x.n <- row.names(x)
    if (is.null(colnames(a))) 
        stop("'a' must have column names.", "\n")
    else a.n <- colnames(a)
    s.x <- dim(x)[1]
    s.a <- dim(a)[2]
    if (s.x != s.a) 
        stop("Different number of species in 'x' and 'a'.", "\n")
    if (any(x.n != a.n)) 
        stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
            "\n")
    com <- dim(a)[1]
    t <- dim(x)[2]
    com.names <- row.names(a)
    sp.names <- row.names(x)
    tr.names <- names(x)
    a[which(is.na(a))] <- 0
    CWM.type <- match.arg(CWM.type)
    is.bin <- function(k) all(k[!is.na(k)] \%in\% c(0, 1))
    bin.var <- rep(NA, t)
    names(bin.var) <- tr.names
    for (i in 1:t) bin.var[i] <- is.bin(x[, i])
    if (!all(bin.var[bin.num])) 
        stop("'bin.num' points to non-binary variables.\n")
    bin.var[bin.num] <- FALSE
    type <- sapply(x, data.class)
    type[type \%in\% c("numeric", "integer")] <- "C"
    type[type == "ordered"] <- "O"
    type[type == "factor"] <- "N"
    type[bin.var] <- "B"
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
        a <- a * biomassValue2
    }
    sum.a <- apply(a, 1, sum)
    a <- a/sum.a
    a <- t(a)
    a <- data.frame(a)
    temp <- list()
    for (i in 1:t) {
        if (type[i] == "C") {
            vec <- numeric(com)
            for (j in 1:com) vec[j] <- weighted.mean(x[, i], 
                a[, j], na.rm = T)
            temp[[i]] <- matrix(vec, com, 1, dimnames = list(com.names, 
                tr.names[i]))
        }
        else {
            x[, i] <- as.factor(x[, i])
            fac <- data.frame()
            which.dom <- rep(NA, com)
            for (k in 1:com) {
                temp2 <- tapply(a[, k], x[, i], sum)
                fac <- rbind(fac, temp2)
                which.dom[k] <- sample(levels(x[, i])[which(fac[k, 
                  ] == max(fac[k, ]))], size = 1)
            }
            colnames(fac) <- paste(tr.names[i], "_", levels(x[, 
                i]), sep = "")
            rownames(fac) <- com.names
            which.dom <- data.frame(which.dom)
            colnames(which.dom) <- tr.names[i]
            rownames(which.dom) <- com.names
            if (CWM.type == "dom") 
                temp[[i]] <- which.dom
            if (CWM.type == "all") 
                temp[[i]] <- fac
        }
    }
    temp <- data.frame(temp)
    return(temp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
