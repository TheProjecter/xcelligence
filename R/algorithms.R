##----------------------------------------------------------------------------##
## Algorithms used in RTCA data analysis
##
## Author: Jitao David Zhang <j.zhang@dkfz.de>
##         Division of Molecular Genome Analysis, DKFZ, Heidelberg, Germany
##
##----------------------------------------------------------------------------##


##--------------------------------------------------##
## transformation methods
##--------------------------------------------------##

###--------------------###
### ratio transform
###--------------------###
ratioTransform <- function(object, time) {
  istime <- nearestTimeIndex(object, time)

  frame <- exprs(object)
  for(i in 1:ncol(frame)) {
    frame[,i] <- frame[,i]/frame[istime,i]
  }
  exprs(object) <- frame
  return(object)
}

###--------------------###
### smoothTransform
###--------------------###
smoothTransform <- function(object,...) {
  time <- timepoints(object)
  frame <- exprs(object)
  for(i in 1:ncol(frame)) {
    if(all(is.na(frame[,i]))) {
      next;
    }
    ss <- smooth.spline(time, frame[,i],...)
    frame[,i] <- ss$y
  }
  exprs(object) <- frame
  return(object)
}

###--------------------###
### interpolationTransform
###--------------------###
.approxRTCAsingle <- function(x, times, newtimes, method="linear") {
  stopifnot(length(x)==length(times))
  if(method %in% c("fmm","periodic","natural","monoH.FC")) {
    res <- spline(times, x, xout=newtimes, method=method)
  } else {
    res <- approx(times,x, xout=newtimes, method=method)
  }
  return(res$y)
}
interpolationTransform <- function(object, interval=0.01, method=c("linear","constant","fmm","periodic","natural", "monoH.FC")) {
  method <- match.arg(method)

  ## extract time points, exprs and new time points
  tps <- timepoints(object)
  readout <- exprs(object)
  itps <- seq(min(tps), max(tps), interval)

  res <- apply(readout, 2, .approxRTCAsingle, times=tps, newtimes=itps, method=method)
  exprs(object) <- res
  timepoints(object) <- itps

  return(object)
}

###--------------------###
### derivativeTransform
###--------------------###
## derivative transform

.Deriv1 <- function(x,y) {
  y.prime <- diff(y) / diff(x)
  x.prime <- x[-length(x)] + diff(x)/2
  list(x = x.prime,
       y = y.prime)
}

derivativeTransform <- function(object) {
  frame <- exprs(object)
  time <- timepoints(object)
  newframe <- frame
  
  for(i in 1:ncol(frame)) {
    y <- frame[,i]
    der <- .Deriv1(time, y)
    newframe[,i] <- c(der$y, der$y[length(der$y)])
  }
  exprs(object) <- newframe
  return(object)
}

###--------------------###
### relative growth rate
###--------------------###
rgrTransform <- function(object, smooth=TRUE) {
  der <- derivativeTransform(object)
  
  deltay <- exprs(der)
  object <- smoothTransform(object)
  y <- exprs(object)
  
  k <- deltay / y
  ## first row replaced by second
  k[1L,] <- k[2L,]
  
  exprs(object) <- k
  if(smooth)
    object <- smoothTransform(object)
  return(object)
}


##--------------------------------------------------##
## statistical methods
##--------------------------------------------------##

## rtcawatershed: watershed algorithm to detect peaks
rtcawatershed <- function(x, thres=quantile(abs(x), 0.9)) {
  x <- abs(x)
  x[x<thres] <- 0
  v <- x >= thres
  vi <- which(v)
  vid <- diff(vi)
  vid <- append(vid, 0) ## peaks at the last time points
  
  r2 <- which(vid != 1) ## breaking point
  if(length(r2)==1) {
    res <- new("rtcaWatershed",
               peakCenter=mean(vi),
               peakInterval=length(vi),
               peakAvgHeight=mean(x[vi]))
  } else {    
    r1 <- append(0, r2[1:(length(r2)-1)])+1

    ui <- mapply(function(x,y) unname(vi[seq(from=x, to=y)]), r1, r2)
    uiMaxLength <- ui[[which.max(sapply(ui, function(x) length(x)))]]

    res <- new("rtcaWatershed",
               peakCenter=mean(uiMaxLength),
               peakInterval=length(uiMaxLength),
               peakAvgHeight=mean(x[uiMaxLength]))
  }
  return(res)
}

