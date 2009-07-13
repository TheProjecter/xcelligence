\name{plotGridEffect}
\Rdversion{1.1}
\alias{plotGridEffect}
\title{
  PLOT GRID EFFECT OF RTCA
}
\description{
  Plot the mean and deviation of rows/columns of a RTCA \emph{E-plate},
  to provide hints of potential row/column effect of the plate
}
\usage{
plotGridEffect(rtca, mode = c("column", "row"), xlab = "time point", ylab = "readout", legend = TRUE, ...)
}
\arguments{
  \item{rtca}{An object of \code{\linkS4class{RTCA}}}
  \item{mode}{character, either \dQuote{column} or \dQuote{row}, to
    choose which effect to depict}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{legend}{logical, whether the legend should be added}
  \item{\dots}{Further parameters passed to \code{\link{plot}} function}
}
\details{
  The error bars depicts the standard deviations
}
\value{
  \code{NULL}, the funciton is called for its side effect
}
\author{Jitao David Zhang}
\examples{
require(RTCA)
  
ofile <- system.file("extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)
plotGridEffect(x)
}
\keyword{ hplot }

