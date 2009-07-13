\name{combineRTCA}
\Rdversion{1.1}
\alias{combineRTCA}
\title{
  Combine a list of RTCA objects
}
\description{
  Combine a list of RTCA objects
}
\usage{
combineRTCA(list)
}
\arguments{
  \item{list}{A list of \code{RTCA} objects}
}
\details{
  The current implementation requires all the objects have exactly the
  same time-points recorded (or at least of same length).

  The combined \code{RTCA} object has an obligatory column in the
  \code{phenoData} \sQuote{Plate} (upper-case!), which matches the names of the
  \code{RTCA} list. When the \code{list} has no names, the
  \sQuote{Plate} field is filled with integer index starting from 1.
}
\note{
  Special attention should be given to the cases where the \code{list}
  parameter partially has names. In this case all items without name
  will be assigned to a \sQuote{Plate} field of empty string
  (\dQuote{}). Therefore it is advised either to assign names to all
  items of the list, or leave them all off.
}
\value{
  A new \code{RTCA} object
}
\author{
  Jitao David Zhang \email{j.zhang@dkfz.de}
}
\examples{
## An artificial example
require(RTCA)
  
ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

xSub1 <- x[,1:3]
xSub2 <- x[,4:ncol(x)]
xComb <- combineRTCA(list(sub1=xSub1, sub2=xSub2))
identical(exprs(x), exprs(xComb))
pData(xComb)$Plate

## in case of nameless list
pData(combineRTCA(list(xSub1, xSub2)))$Plate

## partial names
pData(combineRTCA(list(a=xSub1, xSub2)))$Plate
}
