\name{factor2numeric}
\Rdversion{1.1}
\alias{factor2numeric}
\alias{relevels}
\title{FACTOR UNITILITIES}
\description{
  The functions implement easy interface to certain tasks of factor. See
  datails for explaination
}
\usage{
factor2numeric(x)
relevels(x, refs)
}
\arguments{
  \item{x}{A vector of factor}
  \item{refs}{A vector of character, reference vector to give the orderof levels}
}
\details{
  \code{\link{relevels}} re-arrange the order of levels by the given
  character \code{refs}. Alternatively user could use
  \code{factor(\ldots,levels=refs)} to achieve a similar effect, however
  the \code{\link{relevels}} enables also partial list. The missing
  levels in \code{refs} will be ordered to the last.

  \code{\link{factor2numeric}} converts factor of numerics into their
  numeric form.
}
\value{
  A vector of factor
}
\author{Jitao David Zhang \email{j.zhang@dkfz.de}}
\examples{
## factor2numeric
numFac <- factor(c(3.5, 2.5, 2.5,3.5, 1))
numFac
levels(numFac)

factor2numeric(numFac)
class(factor2numeric(numFac))

## relevels
relevels(numFac, c("3.5", "1", "2.5"))
relevels(numFac, c("3.5", "2.5"))
}
\keyword{ misc }

