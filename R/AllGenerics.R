## RTCAtimeline
setGeneric("actionTrack", function(object) standardGeneric("actionTrack"))
setGeneric("actionTrack<-", function(object, value) standardGeneric("actionTrack<-"))
setGeneric("timeUnit", function(object) standardGeneric("timeUnit"))
setGeneric("timeUnit<-", function(object, value) standardGeneric("timeUnit<-"))
setGeneric("startTime", function(object) standardGeneric("startTime"))
setGeneric("startTime<-", function(object, value) standardGeneric("startTime<-"))
setGeneric("getAction", function(object, time,...) standardGeneric("getAction"))
setGeneric("addAction", function(object, time, action,...) standardGeneric("addAction"))
setGeneric("rmAction", function(object, time, ...) standardGeneric("rmAction"))
setGeneric("orderAction", function(object) standardGeneric("orderAction"))
setGeneric("updateAction", function(object, time, action,...) standardGeneric("updateAction"))
setGeneric("reset", function(object) standardGeneric("reset"))

## RTCAtimeline
setGeneric("timeline", function(object) standardGeneric("timeline"))
setGeneric("timeline<-", function(object, value) standardGeneric("timeline<-"))


## time point
setGeneric("timepoints", function(object) standardGeneric("timepoints"))
setGeneric("timepoints<-", function(object,value) standardGeneric("timepoints<-"))

## plot
setGeneric("plotRTCA", function(x,y,...) standardGeneric("plotRTCA"))
