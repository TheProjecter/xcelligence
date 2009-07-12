\name{rgrTransform}
\Rdversion{1.1}
\alias{rgrTransform}
\title{
  TRANSFORM RTCA DATA INTO RELATIVE GROWTH RATE
}
\description{
  Transform RTCA data into relative growth rate
}
\usage{
rgrTransform(object, smooth)
}
\arguments{
  \item{object}{An object of \code{\linkS4class{RTCA}}}
  \item{smooth}{logical, should the object be smooth transformed after
    the \code{rgrTransform}? Set to \code{TRUE} by default}
}
\details{
  TODO: relative growth rate
}
\value{
  An object of \code{\linkS4class{RTCA}} populated with relative growth
  rate instead of input data
}
\references{
  TODO: reference
}
\author{
  Jitao David Zhang <j.zhang@dkfz.de>
}
\seealso{
  \code{\link{derivativeTransform}} for first derivative. \code{\link{ratioTransform}} for ratio normalization
  adopted by Roche commercial software. \code{\link{smoothTransform}}
  and \code{\link{interpolationTransform}} for other transformation possibilities.
}
\examples{
require(RTCA)
  
ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

xRgr <- rgrTransform(x)
}
\keyword{ models }
