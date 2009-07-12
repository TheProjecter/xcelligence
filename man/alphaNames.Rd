\name{alphaNames}
\Rdversion{1.1}
\alias{alphaNames}
\alias{repairAlphaName}
\alias{alphaNames2Pos}
\alias{rowcol2pos}
\title{
  Auxilliary functions for experiments with microtitre plates
}
\description{
  Functions to manipulate indices or names of microtitre plates
}
\usage{
alphaNames(row = 8, column = 12, order=c("column","row"))
repairAlphaName(x)
alphaNames2Pos(x)
rowcol2pos(row = 1, column=1, plateFormat=c("96","384"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{row}{integer, row index, 1,\ldots,8 for 96-well plates}
  \item{column}{integer, column index, 1,\ldots,12 for 96-well plates}
  \item{x}{character, Well alpha name, in the form of [A-Z][0-9][0-9],
    like 'A01'}
  \item{order}{character, should the alpha names returned in a row-first
    or column-first order?}
  \item{plateFormat}{integer, the microtitre format, either 96 or 384}
  
}
\details{
  \code{alphaNames} returns so-called \emph{alpha well names} in the form of
  [A-H][0-9][0-9] (i.e., A01, C03, D11, H12) for microtitre plates. The
  order of returned alphaNames is controlled by the option \code{order},
  which can be set either as \code{col} or \code{row}

  \code{repairAlphaName} attempts to fix incomplete alpha well names. Now it
  is mainly used to fix well names missing the leading 0 of numeric
  index, like A1.

  \code{alphaName2Pos} returns the row and column number of the given
  alpha well name, in the form of two-column data frame with \emph{row}
  and \emph{col} as colnames.

  \code{rowcol2pos} returns the row-wise position index of given row and
  column index.
}
\value{
  See details
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\examples{
wells <- alphaNames()

repairAlphaName("A1")

alphaNames2Pos(c("A01","B02","C03","H12"))

rowcol2pos(3,1)
}
