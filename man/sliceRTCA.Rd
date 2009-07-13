\name{sliceRTCA}
\Rdversion{1.1}
\alias{sliceRTCA}
\title{
  SLICE RTCA OBJECT WITH TIME
}
\description{
  Subset (slice) RTCA object with starting- and ending-time
}
\usage{
sliceRTCA(x, start, end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of \code{\linkS4class{RTCA}}}
  \item{start}{numeric, start time}
  \item{end}{numeric, end time}
}
\details{
  In case the exact starting- or ending-time is not matched, the nearest
  time point will be used to subset.
}
\value{
  An object of \code{\linkS4class{RTCA}}
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\examples{
  require(RTCA)
  
  ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
  x <- parseRTCA(ofile)

  subx <- sliceRTCA(x, 20, 50)
}