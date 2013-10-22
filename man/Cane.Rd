\name{Cane}
\alias{Cane}
\title{Converts bee intertegular distance (IT; in cm) to body mass (gr)}
\usage{
  Cane(IT)
}
\arguments{
  \item{IT}{A vector of IT measurments \code{IT}}
}
\value{
  output A vector of bee body masses
}
\description{
  The conversion is done as per Cane et al. 1996¿? formula
  and the correlation is >0.96.
}
\examples{
Cane(c(1.2, 2.3, 0.6))
}
\keyword{Bee}

