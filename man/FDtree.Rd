% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/FDtree.R
\name{FDtree}
\alias{FDtree}
\title{Cretae a dendrogram for use in FD analysis}
\usage{
FDtree(S, w = NA, Distance.method = "gower", ord = c("podani", "metric"),
  Cluster.method = c(ward = "ward", single = "single", complete = "complete",
  UPGMA = "average", UPGMC = "centroid", WPGMC = "median", WPGMA = "mcquitty"))
}
\arguments{
\item{S}{matrix or data frame of functional traits. Traits can be numeric, ordered,
or factor. NAs are tolerated.\code{}}

\item{w}{vector listing the weights for the traits in x. Can be missing,
in which case all traits have equal weights.\code{}}

\item{Distance.method}{metric to calculate the species distance matrix. Only Gower is
implemented. \code{}}

\item{ord}{character string specifying the method to be used for ordinal traits
(i.e. ordered). "podani" refers to Eqs. 2a-b of Podani (1999), while "metric"
refers to his Eq. 3. See gowdis for more details.\code{}}

\item{Cluster.method}{Distance method used to produce the tree. UPGMA="average" is
usually giving th ebest results (Podani et al. 2011)\code{}}
}
\value{
an xtree object

The queality of the tree is printed. The quality of the dendogram representation.
clustering  performance is assessed by the correlation with the cophenetic distance
}
\description{
Calculate dendrogram and extract branch lengths for use in FD analysis.
}
\examples{
ex1 <- FDtree(S = dummy$trait, w = NA,
                   Distance.method = "gower", ord = "podani", Cluster.method = "average")
ex1
}

