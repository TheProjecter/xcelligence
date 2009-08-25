setClass("RTCAtimeline",
         representation=representation(actionTrack="data.frame", timeUnit="character", startTime="POSIXct"),
         prototype=prototype(actionTrack=data.frame(time=as.numeric(0), action=I("start")), timeUnit="hour", startTime=Sys.time()),
         validity=function(object) {
           atrack <- actionTrack(object)
           val <- ""
           ## data.frame colnames
           if(all(colnames(atrack)==c("time","action"))) {
             msg <- ""
           } else {
             msg <- "RTCAtimeline must be a two-column data frame, with colnames 'time' and 'action'\n" 

           }
           val <- paste(val, msg, sep="")
           
           ## classes
           if(class(atrack$time) == "numeric" & class(atrack$action) == "character") {
             msg <- ""
           } else {
             msg <- "The two columns must be of class 'numeric' and 'character'\n"
           }
           val <- paste(val, msg, sep="")

           if(val=="")
             val <- TRUE
           
           return(val)
         })

setClass("RTCAhistory",
         representation=representation(history="data.frame"),
         prototype=prototype(history=data.frame(command=I(as.character(NULL)), finished=as.logical(NULL))),
         validity=function(object) {
           history <- object@history
           if(all(colnames(history)==c("command","finished"))) {
             val <- TRUE
           } else {
             val <- paste(val, "RTCAhistory must be a two-column data frame, with colnames 'comamnd' and 'finished'\n", sep="")
           }
           
           if(class(history$command) == "character" & class(history$finished) == "logical") {
             val <- TRUE
           } else {
             val <- paste(val, "The two columns must be of class 'character' and 'logical'\n",sep="")
           }
           
           return(val)
         })
         
setClass("RTCA",
         representation(expID="character",timeline="RTCAtimeline"),
         prototype=list(expID=as.character(NA),timeline=new("RTCAtimeline")),
         contains="ExpressionSet")

setClass("rtcaWatershed",
         representation(peakCenter="numeric",
                        peakInterval="numeric",
                        peakAvgHeight="numeric"))

##----------------------------------------##
## obsolete definition of RTCA
##----------------------------------------##
##setClass("RTCA",
##         representation(normbase="numeric",
##                        smooth="character",
##                        interpolate="character",
##                        timeline="RTCAtimeline"),
##         prototype=list(normbase=-1, smooth="none", interpolate="none",timeline=new("RTCAtimeline")),
##         contains="ExpressionSet")
##
