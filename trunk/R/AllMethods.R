## global variables for show method
.g_wid <- options()$width; .g_decwid <- round(.g_wid*0.60) ## decorator width
.g_decs <- paste(paste(rep("=", .g_decwid),collapse=""),"\n", sep="")

##--------------------##
## setter and getter
##--------------------##
setMethod("actionTrack", "RTCAtimeline", function(object) {
  object@actionTrack
})
setReplaceMethod("actionTrack", c("RTCAtimeline","data.frame"), function(object, value) {
  object@actionTrack <- value
  validObject(object)
  return(object)
})
setMethod("timeUnit", "RTCAtimeline", function(object) {
  object@timeUnit
})
setReplaceMethod("timeUnit", c("RTCAtimeline","character"), function(object, value) {
  object@timeUnit <- value
  return(object)
})
setMethod("startTime", "RTCAtimeline", function(object) {
  object@startTime
})
setReplaceMethod("startTime", c("RTCAtimeline","POSIXct"), function(object, value) {
  object@startTime <- value
  return(object)
})

setMethod("timeline", "RTCA", function(object) {
  object@timeline
})
setReplaceMethod("timeline", "RTCA", function(object,value) {
  object@timeline <- value
  return(object)
})

##--------------------##
## RTCAtimeline: actions
##--------------------##
setMethod("reset", "RTCAtimeline", function(object) {
  nobj <- new("RTCAtimeline")
  timeUnit(nobj) <- timeUnit(object)
  startTime(nobj) <- startTime(nobj)
  return(nobj)
})
setMethod("addAction", c("RTCAtimeline","numeric", "character"), function(object, time, action) {
  aTrack <- actionTrack(object)
  aTrack[nrow(aTrack)+1,] <- c(time, action)

  ## casting: since adding new elements into data.frame will cause time in chr. and action as AsIs
  aTrack$time <- as.numeric(aTrack$time)
  aTrack$action <- as.character(aTrack$action)
  
  actionTrack(object) <- aTrack
  object <- orderAction(object)
  validObject(object)
  return(object)
})
setMethod("addAction", c("RTCA","numeric","character"), function(object, time, action) {
  tl <- timeline(object)
  tl <- addAction(tl, time, action)
  timeline(object) <- tl
  return(object)
})
setMethod("getAction", c("RTCAtimeline", "numeric"), function(object, time, warning=FALSE) {
  aTrack <- actionTrack(object)
  act <- aTrack$action[aTrack$time == time]
  if(length(act)==0) {
    act <- NULL
    if(warning) {
      warning("No action found at the given time point")
    }
  }
  return(act)
})
setMethod("getAction",  c("RTCA", "numeric"), function(object, time, warning=FALSE) {
  getAction(timeline(object), time=time, warning=warning)
})

setMethod("rmAction", c("RTCAtimeline","numeric"), function(object, time, warning=FALSE) {
  aTrack <- actionTrack(object)
  found <- match(time,aTrack$time)
  if(is.na(found)) {
    if(warning) {
      warning("No action found at the given time point")
    }
    return(object)
  } else {
    aTrack <- aTrack[-found,]
    actionTrack(object) <- aTrack
    return(object)
  }
})
setMethod("rmAction", c("RTCA","numeric"), function(object, time, warning=FALSE) {
  tl <- timeline(object)
  tl <- rmAction(tl, time, warning=warning)
  timeline(object) <- tl
  return(object)
})

setMethod("updateAction", c("RTCAtimeline","numeric","character"), function(object, time, action, ifnotfound=c("add","ignore"), warning=FALSE) {
  aTrack <- actionTrack(object)
  found <- match(time, aTrack$time)
  ifnotfound <- match.arg(ifnotfound)
  if(is.na(found)) {
    if(warning) {
      warning("No action found at the given time point")
    }
    if(ifnotfound=="add") {
      object <- addAction(object, time, action)
    } else if (ifnotfound=="ignore") {
    }
  } else {
    aTrack[found,"action"] <- action
    actionTrack(object) <- aTrack
  }

  validObject(object)
  return(object)
})

setMethod("updateAction", c("RTCA","numeric","character"), function(object, time, action, ifnotfound=c("add","ignore"), warning=FALSE) {
  tl <- timeline(object)
  ifnotfound <- match.arg(ifnotfound)
  tl <- updateAction(tl, time, action, ifnotfound=ifnotfound, warning=warning)
  timeline(object) <- tl
  return(object)
})

## internal
setMethod("orderAction", "RTCAtimeline", function(object) {
  aTrack <- actionTrack(object)
  aTrack <- aTrack[order(aTrack$time),]
  actionTrack(object) <- aTrack
  return(object)
})

##----------------------------------------##
## show methods
##----------------------------------------##
setMethod("show","RTCAtimeline", function(object) {
  cat("RTCAtimeline\n")
  cat(.g_decs)
  print(actionTrack(object), row.names=FALSE)
  cat(.g_decs)
  cat("Time unit:", timeUnit(object), "\n")
  cat("RTCA-run start time:", as.character(startTime(object)), "\n")
})

setMethod("show", "RTCA", function(object) {
  show(timeline(object))
  cat("\n")
  show(as(object, "ExpressionSet"))
})

##--------------------##
## time point methods
##--------------------##
setMethod("timepoints", "RTCA", function(object) {
  return(pData(featureData(object))$timepoints)
})
setReplaceMethod("timepoints", "RTCA", function(object,value) {
  featureData(object) <- new("AnnotatedDataFrame",
                             data = data.frame(timepoints=value),
                             varMetadata = data.frame(labelDescription=c("Time points")))
  featureNames(object) <- value
  return(object)
})

## an wrapper to timepoints(x)
setMethod("time", "RTCA", function(x,...) {
  return(timepoints(x))
})

##--------------------##
## plot methods
##--------------------##
setMethod("plotRTCA", "RTCA", function(x,y,...) {
 times <- timepoints(x)
 mat <- exprs(x)[,y,drop=FALSE]
 matplot(times, mat,...)
})


##----------------------------------------##
## transformation
##----------------------------------------##
