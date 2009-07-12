\name{RTCAtimeline-class}
\Rdversion{1.1}
\docType{class}
\alias{RTCAtimeline-class}
\alias{addAction,RTCAtimeline,numeric,character-method}
\alias{actionTrack}
\alias{actionTrack,RTCAtimeline-method}
\alias{actionTrack<-}
\alias{actionTrack<-,RTCAtimeline,data.frame-method}
\alias{addAction}
\alias{addAction,RTCAtimeline,numeric,character-method}
\alias{getAction}
\alias{getAction,RTCAtimeline,numeric-method}
\alias{orderAction}
\alias{orderAction,RTCAtimeline-method}
\alias{reset}
\alias{reset,RTCAtimeline-method}
\alias{rmAction}
\alias{rmAction,RTCAtimeline,numeric-method}
\alias{updateAction}
\alias{updateAction,RTCAtimeline,numeric,character-method}
\alias{timeUnit}
\alias{timeUnit,RTCAtimeline-method}
\alias{timeUnit<-}
\alias{timeUnit<-,RTCAtimeline,character-method}
\alias{startTime}
\alias{startTime,RTCAtimeline-method}
\alias{startTime<-}
\alias{startTime<-,RTCAtimeline,POSIXct-method}

\title{Class "RTCAtimeline"}
\description{Time line of actions performed by the xCelligence device}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("RTCAtimeline")}. However, it is more common to be called
  implicitly by creating an instance of \code{\linkS4class{RTCA}}
  object. 
}
\section{Slots}{
	 \describe{
    \item{\code{actionTrack}:}{Object of class \code{"data.frame"},
      records action track in the form of two-column
      \code{data.frame}. The two columns must have the names
      \sQuote{time} and \sQuote{action}.}
    \item{\code{timeUnit}:}{Object of class \code{"character"},
      recording the unit of time points stored in the \code{actionTrack}
    slot.}
    \item{\code{startTime}:}{Object of class \code{"POSIXct"}, the
      absolute time when the measurement started (at the time point
      \sQuote{0})
    }
  }
}
\section{Methods}{
  \describe{
    \item{addAction}{\code{signature(object = "RTCAtimeline", time =
	"numeric", action = "character")}: ...}
    \item{actionTrack}{\code{signature(object = "RTCAtimeline")}: ... }
    \item{actionTrack<-}{\code{signature(object = "RTCAtimeline", value
	= "data.frame")}: ... }
  }
}
\references{
  \describe{
    \item{1}{\url{http://www.xcelligence.roche.com/}
      introduces \emph{xCelligence} system.}
    \item{2}{\url{http://www.roche-applied-science.com/proddata/gpip/3_8_9_1_1_1.html}
      for brief introduction into RTCA}
    }
}
\author{  Jitao David Zhang \email{j.zhang@dkfz.de} }
\seealso{
  \code{\linkS4class{RTCA}}
}
\examples{
tl <- new("RTCAtimeline")
show(tl)
}
\keyword{classes}
