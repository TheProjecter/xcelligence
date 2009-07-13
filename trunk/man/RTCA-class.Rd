\name{RTCA-class}
\Rdversion{1.1}
\docType{class}
\alias{RTCA-class}
\alias{addAction,RTCA,numeric,character-method}
\alias{derivativeTransform,RTCA-method}
\alias{getAction,RTCA,numeric-method}
\alias{interpolationTransform,RTCA-method}
\alias{plotRTCA,RTCA-method}
\alias{ratioTransform,RTCA-method}
\alias{rgrTransform,RTCA-method}
\alias{rmAction,RTCA,numeric-method}
\alias{RTCAtimeline,RTCA-method}
\alias{RTCAtimeline<-,RTCA-method}
\alias{show,RTCA-method}
\alias{smoothTransform,RTCA,missing-method}
\alias{smoothTransform,RTCA,ANY-method}
\alias{time,RTCA-method}
\alias{expID}
\alias{expID,RTCA-method}
\alias{expID<-}
\alias{expID<-,RTCA-method}
\alias{timeline}
\alias{timeline,RTCA-method}
\alias{timeline<-}
\alias{timeline<-,RTCA-method}
\alias{timepoints}
\alias{timepoints,RTCA-method}
\alias{timepoints<-}
\alias{timepoints<-,RTCA-method}
\alias{updateAction,RTCA,numeric,character-method}

\title{Class "RTCA"}
\description{RTCA object}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RTCA", assayData, phenoData, featureData, experimentData, annotation, exprs, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{expID}:}{Object of class \code{"character"}, experiment ID}
    \item{\code{timeline}:}{Object of class \code{"RTCAtimeline"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{\link[Biobase]{ExpressionSet-class}}, directly.
Class \code{\link[Biobase]{eSet-class}}, by class "ExpressionSet", distance 2.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "ExpressionSet", distance 3.
Class \code{\link[Biobase]{Versioned-class}}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{addAction}{\code{signature(object = "RTCA", time = "numeric", action = "character")}: ... }
    \item{derivativeTransform}{\code{signature(object = "RTCA")}: ... }
    \item{getAction}{\code{signature(object = "RTCA", time = "numeric")}: ... }
    \item{interpolationTransform}{\code{signature(object = "RTCA")}: ... }
    \item{plotRTCA}{\code{signature(x = "RTCA")}: ... }
    \item{ratioTransform}{\code{signature(object = "RTCA")}: ... }
    \item{rgrTransform}{\code{signature(object = "RTCA")}: ... }
    \item{rmAction}{\code{signature(object = "RTCA", time = "numeric")}: ... }
    \item{RTCAtimeline}{\code{signature(object = "RTCA")}: ... }
    \item{RTCAtimeline<-}{\code{signature(object = "RTCA")}: ... }
    \item{show}{\code{signature(object = "RTCA")}: ... }
    \item{smoothTransform}{\code{signature(object = "RTCA", method = "missing")}: ... }
    \item{smoothTransform}{\code{signature(object = "RTCA", method =
	"ANY")}: ...}
    \item{expID}{code{signature(object = "RTCA")}: get Experiment ID}
    \item{expID<-}{code{signature(object = "RTCA", value = "ANY")}: get Experiment ID}
    \item{time}{\code{signature(x = "RTCA")}: ... }
    \item{timeline}{\code{signature(object = "RTCA")}: ... }
    \item{timeline<-}{\code{signature(object = "RTCA")}: ... }
    \item{timepoints}{\code{signature(object = "RTCA")}: ... }
    \item{timepoints<-}{\code{signature(object = "RTCA")}: ... }
    \item{updateAction}{\code{signature(object = "RTCA", time = "numeric", action = "character")}: ... }
  }
}
\references{
  \describe{
    \item{1}{\url{https://www.roche-applied-science.com/sis/xcelligence/index.jsp?id=xcect_000000}
      introduces \emph{xCelligence} system.}
    \item{2}{\url{http://www.roche-applied-science.com/proddata/gpip/3_8_9_1_1_1.html}
      for brief introduction into RTCA}
  }
}
\author{ Jitao David Zhang \email{j.zhang@dkfz.de} }
\examples{
new("RTCA")
}
\keyword{classes}
