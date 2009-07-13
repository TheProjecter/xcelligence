\name{controlView}
\Rdversion{1.1}
\alias{controlView}
\title{
  PLOT CONTROL WELLS IN RTCA DATA
}
\description{
    A convenience function to plot sample wells with control wells on an
  \emph{E-plate} in RTCA system. To use the function the phenoData field
  of the \code{RTCA} object must contain a field named \dQuote{GeneSymbol}.
}
\usage{
controlView(rtca, genesymbol = c("Allstar", "COPB2", "GFP", "mock", "PLK1", "WEE1"), cols, ylim, smooth = FALSE, group = TRUE, ylab = "Normalized cell index", xlab = "Time interval (hour)", drawsd = TRUE, normline = TRUE, ncol = 1, legendpos = "topleft", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rtca}{An object of \code{\linkS4class{RTCA}.}To use the function the phenoData field
  of the \code{RTCA} object must contain a field named \dQuote{GeneSymbol}}
  \item{genesymbol}{character, gene symbols to be plotted. }
  \item{cols}{character, colors used by the provided gene symbols}
  \item{ylim}{y-axis lim}
  \item{smooth}{logical, whether the \code{RTCA} object should be
    smoothed before plotting}
  \item{group}{logical. If \sQuote{group} is set to \code{TRUE}, wells
    with the same \emph{GeneSymbol} will be summarized and plotted. For
    instance, these could be biological replicates. Otherwise each well
    is plotted separatedly}
  \item{ylab}{y axis label}
  \item{xlab}{x axis label}
  \item{drawsd}{logical, should the error bar be drawn to represent
    standard deviation?}
  \item{normline}{logical, should the base-time indicated by a line? See
  \code{\link{ratioTransform}} for the concept of the \emph{base-time}}
  \item{ncol}{integer, legend column number}
  \item{legendpos}{character, legend position}
  \item{\dots}{other parameters passed to the \code{\link{plot}} function}
}
\details{
  The function is often called to draw sample and control in one plot.
}
\value{
  \code{NULL}, the function is called for its side effect
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\seealso{
  \code{\linkS4class{RTCA}}
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

controlView(x, genesymbol=c("mock","COPB2","PLK1"),ylim=c(0,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

