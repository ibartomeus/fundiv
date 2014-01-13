\name{fdisp_w}
\alias{fdisp_w}
\title{Functional Dispersion from FD package with the option to weith by biomass}
\usage{
  fdisp_w(d, a, tol = 1e-07,
    Weigthedby = c("abundance", "biomasCarabids", "biomasBees", "biomassValue"),
    biomassValue = NA)
}
\arguments{
  \item{same}{as dbFD\code{}}

  \item{Weigthedby}{character string indicating if should
  be weighted by `abundance` or `biomassValue`. If
  biomassValue is in "length" units for Carabids or bees,
  use options `biomasCarabids` or `biomasBees` to
  automatically convert to mass.\code{}}

  \item{biomassValue}{numerical vector with body weigh (or
  length) values for each species in the same order as
  species are provided.  \code{}}
}
\value{
  same as dbFD\code{}

  FDis vector listing the FDis of each community weighted
  by abundance or biomass
}
\description{
  fdisp measures the functional dispersion (FDis) of a set
  of communities, as described by Laliberté and Legendre
  (2010). fdisp_w also allows to weight FDis by biomass
}
\examples{
dummy.dist <- gowdis(dummy$trait)
ex1 <- fdisp_w(dummy.dist, dummy$abun, Weigthedby = "biomassValue",
biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))
ex1
}

