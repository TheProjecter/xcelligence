\name{plateView}
\Rdversion{1.1}
\alias{plateView}
\title{
  PLATE VIEW OF RTCA DATA
}
\description{
  Plots a \emph{E-plate} in RTCA assays in one plot to convey an
  overview of the plate
}
\usage{
plateView(rtca, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtca}{An object of \code{\linkS4class{RTCA}}}
  \item{ylim}{ylab lim}
  \item{\dots}{Other parameters passed to \code{\link{plot}} function}
}
\details{
  For now the function only supports the visualization of a 96-well \emph{E-plate}.
}
\value{\code{NULL}, the function is called for the side effect}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\seealso{
  \code{\link{plotRTCA}}
}
\examples{
require(RTCA)
  
ofile <- system.file("extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

## Not run automatically, because of 'margin too large'
## plateView(x)
}

\keyword{hplot}
