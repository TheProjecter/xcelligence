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
  A csv-like format file can be exported from the RTCA device, which can
  be fed into this function to set up an instance of
  \code{\linkS4class{RTCA}} object.

  In the \emph{/extdata/} directory of the package, such a file is
  provided as an example. The first line contains the experiment ID,
  which is followed by a matrix of recorded data in the tabular
  form. The first and second column records the time-interval in the
  unit of hour and hour-minute-second format respectively. The rest
  columns then record the read-out (\sQuote{Cell-Index}, or \sQuote{CI})
  of the device, with each well a role.

  \code{phenoData} allows user to annotate the wells. \code{skipWell}
  allows to skip wells in case, for example, they are known to be
  contaminated.

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
  
ofile <- system.file("extdata/testOutput.csv", package="RTCA")
pfile <- system.file("extdata/testOutputPhenoData.csv", package="RTCA")

pData <- read.csv(pfile, sep="\t", row.names="Well")
metaData <- data.frame(labelDescription=c(
"Rack number",
"siRNA catalogue number",
"siRNA gene symbol",
"siRNA EntrezGene ID",
"siRNA targeting accession"
))

phData <- new("AnnotatedDataFrame", data=pData, varMetadata=metaData)
x <- parseRTCA(ofile, phenoData=phData)

x
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ file }

