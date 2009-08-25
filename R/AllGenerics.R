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

## RTCA
setGeneric("expID", function(object) standardGeneric("expID"))
setGeneric("expID<-", function(object,value) standardGeneric("expID<-"))
setGeneric("timeline", function(object) standardGeneric("timeline"))
setGeneric("timeline<-", function(object, value) standardGeneric("timeline<-"))


## time point
setGeneric("timepoints", function(object) standardGeneric("timepoints"))
setGeneric("timepoints<-", function(object,value) standardGeneric("timepoints<-"))

## watershed
setGeneric("peakCenter", function(object) standardGeneric("peakCenter"))
setGeneric("peakCenter<-", function(object,value) standardGeneric("peakCenter<-"))
setGeneric("peakInterval", function(object) standardGeneric("peakInterval"))
setGeneric("peakInterval<-", function(object,value) standardGeneric("peakInterval<-"))
setGeneric("peakAvgHeight", function(object) standardGeneric("peakAvgHeight"))
setGeneric("peakAvgHeight<-", function(object,value) standardGeneric("peakAvgHeight<-"))

## plot
setGeneric("growthHeatmap", function(object,...) standardGeneric("growthHeatmap"))
