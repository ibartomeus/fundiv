\name{length_to_mass}
\alias{Cane}
\alias{Jelaska}
\alias{length_to_mass}
\title{Converts length measures to body mass for bees and Carabids}
\usage{
  length_to_mass(x, taxa = c("Bees", "Carabids"))

  Cane(x)

  Jelaska(x)
}
\arguments{
  \item{x}{A vector of bee intertegular spans (IT)
  measurments in cm or carabid beetles body length
  (mm)\code{}}

  \item{taxa}{string character indicating if length measure
  is from "Bees" or "Carabids"}
}
\value{
  A vector of bee body masses is returned in (gr) for bees
  and (mg) for Carabids.
}
\description{
  For Bees, the conversion is done as per Cane et al. 1987
  (Estimation of bee size using intertegular span
  (Apoidea). Journal of the Kansas Entomological Society
  60:145–147). The correlation between IT and body mass is
  0.96. For Carabids the conversion is done from Jelaska et
  al. 2011 (Carabid beetle diversity and mean individual
  biomass in beech forests of various ages. ZooKeys
  100:393-405) which in turn is based on Szyszko J (1983)
  book.
}
\examples{
length_to_mass(c(1.2, 2.3, 0.6), "Bees")
Jelaska(c(1.2, 2.3, 0.6))
}
\keyword{Bee,}
\keyword{Carabids}

