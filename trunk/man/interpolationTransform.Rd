\name{interpolationTransform}
\Rdversion{1.1}
\alias{interpolationTransform}
\title{
  TRANSFORM RTCA DATA WITH INTERPOLATION
}
\description{
  Interpolate RTCA data
}
\usage{
interpolationTransform(object, interval=0.01, method=c("linear","constant","fmm","periodic","natural", "monoH.FC"))
}
\arguments{
  \item{object}{An \code{RTCA} object}
  \item{\dots}{other parameters, \code{interval} and \code{method} are
    implemented, see below}
  \item{interval}{numeric, the interval between interpolated points, set
    to 0.01 by default}
  \item{method}{character, specifying the method for interpolation,
    \dQuote{linear} by default (for linear interpolation). Allowed options are: \dQuote{linear} and
    \dQuote{constant} for \code{approx} interpolation, and \dQuote{fmm},
    \dQuote{periodic}, \dQuote{natural} and \dQuote{monoH.FC} for cubic
    spline interpolation}

}
\details{
  Since most RTCA experiements record the experiments in the irregular
  time-series, sometimes however it is desired to have regular
  intervals. \code{interpolationTransform} interpolate between data
  points to estimate results of regular intervals.

  Two classes of interpolations are supported by now: linear (using
  \code{\link[stats]{approx}}) and cubic spline
  (\code{\link[stats]{spline}}) interpolation. By default linear
  interpolation is used.
}
\value{
  An interpolated object of \code{\linkS4class{RTCA}}.
}
\author{Jitao David Zhang \email{j.zhang@dkfz.de}}

\seealso{
\code{\link{rgrTransform}} stands for \emph{relative growth rate
  transformation}, \code{\link{ratioTransform}} for ratio normalization
  adopted by Roche commercial software. \code{\link{smoothTransform}} to
  smooth the RTCA readout.
}
\examples{
require(RTCA)
  
ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

xInter <- interpolationTransform(x)
}
\keyword{ ts }

