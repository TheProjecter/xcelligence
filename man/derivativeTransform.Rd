\name{derivativeTransform}
\Rdversion{1.1}
\alias{derivativeTransform}
\title{
  DERIVATIVE TRANSFORM OF RTCA OBJECT
}
\description{
  Derivative transform of RTCA object, returning the change rate of cell impedance
}
\usage{
derivativeTransform(object)
}
\arguments{
  \item{object}{An object of \code{\linkS4class{RTCA}}}
}
\details{
  The first derivative of the cell impedance curve measured by
  RTCA. The derivative of the last time point is estimated by that of
  the next to last point.
}
\value{
  An \code{\linkS4class{RTCA}} object populated with derivative values
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\seealso{
  \code{\link{smoothTransform}} and \code{\link{interpolationTransform}}
  for smoothing and interpolating the RTCA
  data. \code{\link{rgrTransform}} calculates relative growth rate,
  which calls \code{derivativeTransform}.
}
\examples{
require(RTCA)
  
ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

xDeriv <- derivativeTransform(x)
}
