\name{Eve}
\alias{Eve}
\title{Calculate evenness indexes including popular ones and renyi profiles.}
\usage{
Eve(A, scales = c(1))
}
\arguments{
  \item{A}{matrix containing the abundances of the species
  in S (or presence-absence, i.e. 0 or 1). Rows are sites
  and species are columns. NA not tolerated. The number of
  species (columns) in A must match the number of species
  (rows) in S. In addition, the species labels in A and S
  must be identical and in the same order.\code{}}

  \item{scales}{vector containing the values of a (hill
  numbers) to be calculated. The default is 1, giving you
  E1,0 index proposed by Riccota and Avena 2003. A more
  inclusive exploration would be c(0,0.25,0.5,1,2,4,8,Inf)}
}
\value{
Ea0 returns a data.frame of class renyi with the selected
indices.

EinvD returns a vector with EinvD.

EJ returns a a vector with Pielou's index.

Evar returns a vector with Evar.
}
\description{
Calculate eveness indexes for a set of communities using
scripts and the vegan package and returns a single
dataframe. After reading a tone, Pielou's J is included for
historical reasons. Evar and InvSimpson are recomended in X
and I Riccota and Avena 2003 recomend to report the
spectrum of hill numbers.
}
\examples{
Ea0(A = dummy$abun)
Ea0(A = dummy$abun, scales = c(0.25,0.5,1,2,4,8,Inf))
EinvD(A = dummy$abun)
EJ(A = dummy$abun)
Evar(A = dummy$abun)
#calculate all of them
eve1 <- Eve(A = dummy$abun)
eve2 <- Eve(A = dummy$abun, scales = c(0.25,0.5,1,2,4,8,Inf))
pairs(eve1)
pairs(eve2)
}

