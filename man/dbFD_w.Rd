\name{dbFD_w}
\alias{dbFD_w}
\title{Distance-Based Functional Diversity Indices from FD package with the option to weith by biomass}
\usage{
dbFD_w(x, a, w, w.abun = TRUE, stand.x = TRUE, ord = c("podani",
  "metric"), asym.bin = NULL, corr = c("sqrt", "cailliez", "lingoes",
  "none"), calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
  scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
  km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100,
  km.crit = c("calinski", "ssi"), calc.CWM = TRUE, CWM.type = c("dom",
  "all"), calc.FDiv = TRUE, dist.bin = 2, print.pco = FALSE,
  messages = TRUE, Weigthedby = c("abundance", "biomasCarabids",
  "biomasBees", "biomassValue"), biomassValue = NA)
}
\arguments{
  \item{same}{as dbFD\code{}}

  \item{Weigthedby}{character string indicating if should
  be weighted by `abundance` or `biomassValue`. When
  biomassValue is in "length" units for Carabids or Bees,
  use options `bi09:43:19omasCarabids` or `biomasBees` to
  automatically convert it to mass.\code{}}

  \item{biomassValue}{numerical vector with body weigh (or
  length) values for each species in the same order as
  species are provided. It can also be a matrix or data
  frame with one mass value for each community and species
  (both communities and species arranged like in a.
  \code{}}
}
\value{
same as dbFD\code{}

FDis vector listing the FDis of each community weighted by
abundance or biomass

RaoQ vector listing the RaoQ of each community weighted by
abundance or biomass

FEve vector listing the FEve of each community weighted by
abundance or biomass

FDiv vector listing the FDiv of each community weighted by
abundance or biomass
}
\description{
dbFD implements a flexible distance-based framework to
compute multidimensional functional diversity (FD) indices.
dbFD returns the three FD indices of Villéger et al.
(2008): functional richness (FRic), functional evenness
(FEve), and functional divergence (FDiv), as well
functional dispersion (FDis; Laliberté and Legendre 2010),
Rao's quadratic entropy (Q) (Botta-Dukát 2005), a
posteriori functional group richness (FGR) (Petchey and
Gaston 2006), and the community-level weighted means of
trait values (CWM; e.g. Lavorel et al. 2008). Some of these
FD indices consider species abundances. dbFD includes
several options for flexibility. dbFD_w also allows to
weight all metrics by biomass
}
\examples{
library(FD)
ex1 <- dbFD_w(x = dummy$trait, a = dummy$abun, Weigthedby = "biomassValue",
biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))
ex1
}

