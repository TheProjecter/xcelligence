\name{smoothTransform}
\Rdversion{1.1}
\alias{smoothTransform}
\title{
  SMOOTH TRANSFORM OF RTCA OBJECT
}
\description{
  Smoothing the RTCA cell impedance measurement
}
\usage{
smoothTransform(object, ...)
}
\arguments{
  \item{object}{An object of \code{\linkS4class{RTCA}}}
  \item{\dots}{Parameters passed to \code{\link[stats]{smooth.spline}}}
}
\details{
  \code{smoothTransform} smooths the RTCA cell impedance measurement by
  calling the function \code{\link[stats]{smooth.spline}}. This feature
  can be useful for visualiation purposes and in conjuction with other transformations.
}
\value{
   An \code{\linkS4class{RTCA}} object populated with smoothed values
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\note{
  \code{\link{ratioTransform}} performs ratio transformation recommended
  by the machine provider. \code{\link{interpolationTransform}}
  for interpolating the RTCA
  data. \code{\link{derivativeTransform}} returns cell impedance change
  rates and \code{\link{rgrTransform}} calculates relative growth rate.

}
\examples{
require(RTCA)
  
ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

xSmooth <- smoothTransform(x)
}