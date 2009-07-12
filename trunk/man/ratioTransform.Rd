\name{ratioTransform}
\Rdversion{1.1}
\alias{ratioTransform}
\title{
  RATIO TRANSFORMATION OF RTCA DATA
}
\description{
  Performs ratio transformation (normalisation) of RTCA data, as
  recommended by the producer Roche.
}
\usage{
ratioTransform(object, time)
}
\arguments{
  \item{object}{An object of \code{\linkS4class{RTCA}}}
  \item{time}{numeric, the time point used to normalize the whole series
  of data}
}
\details{
  The \emph{xCelligence} software provided by Roche performs ratio
  transform implicitly by dividing the time-series impedance measurement
  by the value of a selected time point (so-called 'base-time'), for instance 5 hours after
  compound transfection, in each cell. The aim of this transformation was to
  scale (normalize) the data of different wells, since the normalized
  values of all wells are uniformly 1 at the base-time.

  However, this method is vulnerable to arbitrary selection of the time
  point chosen to normalize. It may be helpful to try several base-time
  values before comparing normalized results.

  See \code{\link{derivativeTransform}} and \code{\link{rgrTransform}}
  for other normalization (scaling) possibilities.
}
\value{
  An object of \code{\linkS4class{RTCA}}, populated with normalized
  value. The normalized values of all wells are uniformly 1 at the base-time.
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\seealso{
   \code{\link{smoothTransform}} and \code{\link{interpolationTransform}}
  for smoothing and interpolating the RTCA
  data. \code{\link{rgrTransform}} calculates relative growth rate,
  \code{derivativeTransform} calculates derivatitve. The later two
  methods are not sensative to the selection of base-time point.
}
\examples{
require(RTCA)
  
ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

xNorm <- ratioTransform(x, 35)
}
