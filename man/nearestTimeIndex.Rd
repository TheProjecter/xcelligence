\name{nearestTimeIndex}
\Rdversion{1.1}
\alias{nearestTimeIndex}
\title{
  GET INDEX FOR NEAREST TIME
}
\description{
  Get index for the nearest time point to the given one. Called internally in many
  time-point related functions.
}
\usage{
nearestTimeIndex(rtca, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtca}{An object of \code{\linkS4class{RTCA}}}
  \item{time}{numeric, a time point}
}
\details{
  The function finds the time point with minimum absolute difference to
  the given time and returns its index.
}
\value{
  An integer, the index of the nearest time point
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\seealso{
  \code{\link{timepoints}} to return all time points of an
  \code{\linkS4class{RTCA}} object.
}
\examples{
require(RTCA)
  
ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

x
xIndex <- nearestTimeIndex(x, 25)
timepoints(x)[xIndex]
}
