\name{parseRTCA}
\Rdversion{1.1}
\alias{parseRTCA}
\title{
  Parse RTCA output file
}
\description{
  The function parses RTCA output file into \code{RTCA} object
}
\usage{
parseRTCA(file, dec = ".", phenoData, skipWell, ...)
}
\arguments{
  \item{file}{character, name of the RTCA output file}
  \item{dec}{decimal sign of the file}
  \item{phenoData}{phenoData}
  \item{skipWell}{character, well(s) to skip}
  \item{\dots}{other parameters passed to \code{\link{read.table}}}
}
\details{
  TODO: To be implemented
}
\value{
  An object of \code{RTCA-class}
}
\references{
\url{http://www.roche-applied-science.com/proddata/gpip/3_8_9_1_1_1.html}
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\examples{
require(RTCA)
  
ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

x
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ file }

