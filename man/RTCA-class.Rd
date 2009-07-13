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
Objects can be created by calls of the form \code{new("RTCA", assayData,
  phenoData, featureData, experimentData, annotation, exprs,
  ...)}. However, it is more common to be constructed by
\code{\link{parseRTCA}} function by reading in RTCA output data directly.
}
\section{Slots}{
  \describe{
    \item{\code{expID}:}{Object of class \code{"character"}, experiment ID}
    \item{\code{timeline}:}{Object of class \code{"RTCAtimeline"},
      recording action track along the time line}
    \item{\code{assayData}:}{Object of class \code{"AssayData"}, assay
      data inherited from \code{ExpressionSet-class}}
    \item{\code{phenoData}:}{Object of class
      \code{"AnnotatedDataFrame"}, pheno data of the assay, annotating
      the wells}
    \item{\code{featureData}:}{Object of class
      \code{"AnnotatedDataFrame"}, feature data of the assay, preserved
      for time-line recording by the package}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}, idle}
    \item{\code{annotation}:}{Object of class \code{"character"}, idle}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"},idle}
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
    \item{addAction}{\code{signature(object = "RTCA", time = "numeric",
	action = "character")}: add action at the specified time, passed
    to the \code{RTCAtimeline} slot}
    \item{getAction}{\code{signature(object = "RTCA", time =
	"numeric")}: get action at the specified time, passed
    to the \code{RTCAtimeline} slot}
    \item{plotRTCA}{\code{signature(x = "RTCA")}: plot RTCA }
    \item{rmAction}{\code{signature(object = "RTCA", time = "numeric")}:
    remove action at the specified time, passed
    to the \code{RTCAtimeline} slot}
    \item{show}{\code{signature(object = "RTCA")}: print method }
    \item{expID}{code{signature(object = "RTCA")}: get Experiment ID}
    \item{expID<-}{code{signature(object = "RTCA", value = "ANY")}: set Experiment ID}
    \item{time}{\code{signature(x = "RTCA")}: deprecated }
    \item{timeline}{\code{signature(object = "RTCA")}: get the \code{RTCAtimeline} slot }
    \item{timeline<-}{\code{signature(object = "RTCA")}: assign the
      \code{RTCAtimeline} slot }
    \item{timepoints}{\code{signature(object = "RTCA")}: get the
      recording time points in a vector }
    \item{timepoints<-}{\code{signature(object = "RTCA")}: assign the
      recording time points }
    \item{updateAction}{\code{signature(object = "RTCA", time =
	"numeric", action = "character")}: update the action at the specified time, passed
    to the \code{RTCAtimeline} slot }
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
new("RTCA", expID="testExp01")
}
\keyword{classes}
