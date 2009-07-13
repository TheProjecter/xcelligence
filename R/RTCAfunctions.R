##----------------------------------------##
## microtitre plate functions
##----------------------------------------##
alphaNames <- function(row=8, column=12, order=c("column","row")) {
  order <- match.arg(order)
  rowalpha <- LETTERS[1:row]
  columnalpha <- seq(1,to=column)
  fs <- "%s%02d"
  if(order == "column") {
    alphas <- mapply(function(x,y) {sprintf(fs, x, y)} , rep(rowalpha, column), rep(columnalpha, each=row))
  } else if (order == "row") {
    alphas <- mapply(function(x,y) {sprintf(fs, x, y)} , rep(rowalpha, each=column), rep(columnalpha, row))
  }
  
  return(unname(alphas))
}

## repairAlphaName: vector friendly
repairAlphaName <- function(x) {
  x <- gsub("\\s","",x)

  xnchar <- nchar(x)
  is2char <- xnchar == 2
  x[is2char] <- paste(substr(x[is2char],1,1), substr(x[is2char],2,2),sep="0")

  isvalid <- x == grep("[A-Z][0-9][0-9]", x, value=TRUE)
  if(all(isvalid)) {
    return(x)
  } else {
    stop("The following names are not repairable!\n", x[!isvalid])
  }
}

## convert alpha name to position index (integer)
alphaNames2Pos <- function(x) {
  nchars <- sapply(x, nchar)
  stopifnot(all(nchars==3))
  rowl <- sapply(x, function(x) grep(substring(x,1,1),LETTERS))
  coll <- sapply(x, function(x) grep(substring(x, 2,3), sprintf("%02d",1:50)))
  df <- data.frame(row=rowl, column=coll)
  rownames(df) <- x
  return(df)
}

rowcol2pos <- function(row=1, column=1, plateFormat=c("96", "384")) {
  plateFormat <- match.arg(plateFormat)
  if(plateFormat=="96") {
    nRow <- 12
  } else if (plateFormat=="384") {
    nRow <- 24
  }
  posr <- (row-1)*nRow
  posc <- column
  pos <- posr+posc
  return(pos)
}

##----------------------------------------##
## factor functions
##----------------------------------------##
factor2numeric <- function(x) {
  n <- as.numeric(levels(x))[x]
  return(n)
}

## relevel factors. Alternatively user could use factor(..., levels=refs). However,
## relevels also receive partial list.
relevels <- function (x, refs) 
{
    if (!all(refs %in% levels(x))) {
        missing <- which(!(refs %in% levels(x)))
        stop("The following levels are not found in x:\n", paste(refs[missing], 
            sep = ","))
    }
    refs <- rev(refs)
    for (i in refs) {
        x <- relevel(x, ref = i)
    }
    return(x)
}


##----------------------------------------##
## file I/O
##----------------------------------------##
parseRTCA <- function(file,dec=".", phenoData, skipWell,...) {
  scans <- scan(file, what=character(0), sep="\n")
  ## experimentID
  expIdIndex <- grep("Experiment ID",scans)
  expId <- gsub(".*ID:\\s*([[:alnum:]]*)[[:space:]]*", "\\1", scans[expIdIndex], extended=TRUE)
  
  skipnum <- grep("^Time", scans)-1
  dt <- read.table(file, skip=skipnum, sep="\t",head=TRUE,dec=dec,...)
  dt <- dt[-1,] ## 0 is doubled
  dt <- dt[,-2] ## remove time interval
  rownames(dt) <- dt[,1]

  tintervals <- dt[,1]
  if(is(tintervals, "factor")) {
    tintervals <- factor2numeric(tintervals)
  }
  
  dt <- dt[,-1]
  colnames(dt) <- gsub("Y\\.\\.","",colnames(dt))
  colnames(dt) <- gsub("\\.$","",colnames(dt))
  twochar <- sapply(colnames(dt), nchar) == 2
  colnames(dt)[twochar] <- repairAlphaName(colnames(dt)[twochar])

  stopifnot(length(tintervals) == nrow(dt))

  ## abnormal
    if(!missing(skipWell)) {
      abgrep <- grep(skipWell, colnames(dt))
      if(length(abgrep) > 0) {
        for(i in seq(along=abgrep)) {
          dt[,abgrep[i]] <- rep(NA, nrow(dt))
        }
      }
    }

  x <- new("RTCA", expID=expId)
  exprs(x) <- as.matrix(dt)
  timepoints(x) <- tintervals
  if(missing(phenoData)) {
    phenoData <- new("AnnotatedDataFrame", data=data.frame(Well=alphaNames(), GeneSymbol=""))
  }
  phenoData(x) <- phenoData
  return(x)
}

##----------------------------------------##
## RTCA object manipulation
##----------------------------------------##
combineRTCA <- function(list) {
  exprss <- lapply(list, exprs)
  pdatas <- lapply(list, pData)
  pnames <- names(list)
  if(is.null(pnames)) {
    pnames <- seq(along=list)
  } 

  for(i in seq(along=pdatas)) {
    pdatas[[i]]$Plate <- pnames[i]
  }
  
  newexprs <- do.call(cbind, exprss)
  newpdata <- do.call(rbind, pdatas)
  newpdata$Plate <- factor(newpdata$Plate)
  
  newobj <- list[[1]]
  exprs(newobj) <- newexprs
  pData(newobj) <- newpdata
  
  return(newobj)
}


nearestTimeIndex <- function(rtca, time) {
  alltime <- timepoints(rtca)
  nearest <- which.min(abs(time -alltime))
  return(nearest)
}

## cut by time
sliceRTCA <- function(x, start, end) {
  tst <- nearestTimeIndex(x, start)
  tend <- nearestTimeIndex(x, end)
  return(x[tst:tend,])
}


##----------------------------------------##
## Transformation
##----------------------------------------##

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
  tps <- timepoints(object)
  readout <- exprs(object)
  itps <- seq(min(tps), max(tps), interval)
  imatch <- match(tps, itps)
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


##----------------------------------------##
## Visualization
##----------------------------------------##
## special function to determine grid effect
.gridEffect <- function(rtca, mode=c("row","col")) {
  mode <- match.arg(mode)
  frame <- exprs(rtca)
  pos <- alphaNames2Pos(colnames(frame))
  if(mode=="row") {
    posl <- as.factor(pos$row)
  } else {
    posl <- as.factor(pos$col)
  }
  gridmean <- gridsd <- matrix(as.numeric(NA), ncol=nlevels(posl), nrow=nrow(frame))
  for(i in 1:nlevels(posl)) {
    ma <- posl == levels(posl)[i]
    gridmean[,i] <- apply(frame[,ma], 1, mean)
    gridsd[,i] <- apply(frame[,ma], 1, sd)
  }
  return(list(mean=gridmean, sd=gridsd))
}

plotGridEffect <- function(rtca, mode=c("col","row"), xlab="time point", ylab="readout",legend=TRUE,...) {
  mode <- match.arg(mode)
  grid <- .gridEffect(rtca, mode)
  mean <- grid[["mean"]]
  sd <- grid[["sd"]]
  require(RColorBrewer)
  cols <- brewer.pal(ncol(mean),"Set3")
  plot(1:nrow(mean), rep(max(mean), nrow(mean)), type="n", ylim=c(min(mean,na.rm=T), max(mean,na.rm=T)*1.2),
       ylab=ylab, xlab=xlab,...)
  for(i in 1:ncol(mean)) {
    for(j in 1:nrow(mean)) {
      segments(j, mean[j,i]-sd[j,i],j, mean[j,i]+sd[j,i], col=cols[i], lwd=0.5)
    }
  }
  for(i in 1:ncol(mean)) {
    lines(1:nrow(mean), mean[,i], col=cols[i], lwd=4)
  }
#  browser()
  if(legend) {
    legend("topright",
           legend=1:ncol(mean), lty=1, lwd=4, col=cols, ncol=4, bty="n")
  }
}

controlview <- function(rtca,
                        genesymbol=c("Allstar","COPB2","GFP","mock", "PLK1","WEE1"),
                        controlcols,
                        ylim,smooth=FALSE, group=TRUE, ylab="Normalized cell index", xlab="Time interval (hour)",  drawsd=TRUE, normline=TRUE, ncol=1, legendpos="topleft",...) {
  timeint <- timepoints(rtca)
  frame <- exprs(rtca)
  pdata <- pData(rtca)
  if(missing(ylim))
    ylim <- quantile(frame, c(0.05, 0.95))
  plot(timeint, frame[,1], type="n", ylim=ylim, xlab=xlab, ylab=ylab,...)
  genesymbols <- relevels(factor(genesymbol), genesymbol)
  if(missing(controlcols))
    controlcols <- brewer.pal(nlevels(genesymbols),"Set2")
  for (i in 1:nlevels(genesymbols)) {
    gs <- levels(genesymbols)[i]
    gsindex <- grep(paste("^",gs,"$",sep=""),pdata$GeneSymbol)
    if(group) {
      gscols <- frame[,gsindex, drop=FALSE]
      gscols <- gscols[,!apply(gscols, 2, function(x) all(is.na(x))), drop=FALSE]
      groupmean <- rowMeans(gscols)
      groupsd <- apply(gscols,1,sd,na.rm=T)/sqrt(length(gsindex))*1.96
      lines(timeint, groupmean, col=controlcols[i], lwd=3)
      if(drawsd) {
        segments(timeint, groupmean+groupsd, timeint, groupmean-groupsd, col=controlcols[i], lwd=1)
        segments(timeint-0.25, groupmean-groupsd, timeint+0.25, groupmean-groupsd, col=controlcols[i])
        segments(timeint-0.25, groupmean+groupsd, timeint+0.25, groupmean+groupsd, col=controlcols[i])
      }
    } else {
      for (j in gsindex) {
        y <- frame[,j]
        if(smooth) y <- smooth(y)
        lines(timeint, y, col=controlcols[i], lwd=3)
      }
    }
  }
  grid()
  legend(legendpos, legend=levels(genesymbols), col=controlcols, lwd=4, cex=1.1, bty="n", ncol=ncol)


  if(normline)
    abline(h=0, col="red", lwd=2, lty=2)
  normpoint <- apply(frame, 1,function(x) all(x==1,na.rm=TRUE))
  if(any(normpoint) && normline) {
    abline(v=timeint[which(normpoint)[1]], lwd=3, lty=5)
  }
}

## display the curve of a plate in one figure
plateView <- function(rtca,ylim,...) {
  wellno <- ncol(rtca)
  stopifnot(wellno == 96) # only supports 96-well for now
  layout(matrix(1:96, nrow=8, ncol=12, byrow=TRUE))
  opar <- par(mgp=c(0,0,0), mar=c(0,1,1.5,0))
  express <- exprs(rtca)
  if(missing(ylim)) {
    ylim <- c(quantile(express,0.02), quantile(express,0.98))
  }
  for (i in 1:96) {
    plot(timepoints(rtca),express[,i], type="l", xlab="", ylab="", axes=FALSE,ylim=ylim,...)
    title(paste(pData(rtca)$Well[i], pData(rtca)[,2][i]))
    abline(h=0, lty=2, col="darkgrey")
  }
  par(opar)
}


##----------------------------------------##
## gene specific functions
##----------------------------------------##


##----------------------------------------##
## defunted functions
##----------------------------------------##

##reNormalize <- function(rtca, time) {
##  rtcatp <- timepoints(rtca)
##  
##  istime <- match(time, rtcatp)
##
##  if(is.na(istime)) {
##    close <- which.min(abs(as.numeric(time) - rtcatp))
##    istime <- close
##  }
##
##  frame <- exprs(rtca)
##  for(i in 1:ncol(frame)) {
##    frame[,i] <- frame[,i]/frame[istime,i]
##  }
##  exprs(rtca) <- frame
##  return(rtca)
##}

##smoothRTCA <- function(rtca) {
##  time <- timepoints(rtca)
##  frame <- exprs(rtca)
##  for(i in 1:ncol(frame)) {
##    if(all(is.na(frame[,i]))) {
##      next;
##    }
##    ss <- smooth.spline(time, frame[,i])
##    frame[,i] <- ss$y
##  }
##  exprs(rtca) <- frame
##  return(rtca)
##}

##interpolateRTCA <- function(x, interval=0.01, method="linear") {
##  tps <- timepoints(x)
##  readout <- exprs(x)
##  itps <- seq(min(tps), max(tps), interval)
##  imatch <- match(tps, itps)
##  res <- apply(readout, 2, approxRTCAsingle, times=tps, newtimes=itps, method=method)
##  exprs(x) <- res
##  timepoints(x) <- itps
##
##  return(x)
##}

##relGrowthRate <- function(rtca, smooth=FALSE) {
##  ## derivative
##  der <- derivative(rtca)
##
##  deltay <- exprs(der)
##  rtca <- smoothRTCA(rtca)
##  y <- exprs(rtca)
##  
##  k <- deltay / y
##  ## first row replaced by second
##  k[1L,] <- k[2L,]
##  
##  res <- rtca
##  exprs(res) <- k
##  if(smooth)
##    res <- smoothRTCA(res)
##  return(res)
##}


##derivative <- function(rtca) {
##  frame <- exprs(rtca)
##  time <- timepoints(rtca)
##  newframe <- frame
##  
##  for(i in 1:ncol(frame)) {
##    y <- frame[,i]
##    der <- Deriv1(time, y)
##    newframe[,i] <- c(der$y, der$y[length(der$y)])
##  }
##  exprs(rtca) <- newframe
##  return(rtca)
##}

##unionDF <- function(x,y) {
##  .Deprecated("unionDF function deprecated")
##  xcols <- sort(colnames(x))
##  ycols <- sort(colnames(y))
##  if(!all(xcols==ycols)) stop()
##  res <- rbind(x[,xcols], y[,ycols])
##
##  ## x decides the order
##  reord <- match(colnames(x), xcols)
##  return(res[,reord])
##}

##gene2exprs <- function(rtca, gene) {
##  subset <- exprs(rtca)[,pData(rtca)$GeneSymbol==gene]
##  return(subset)
##}
##
##timettest <- function(rtca, gene1, gene2, ...) {
##  g1 <- gene2exprs(rtca, gene1)
##  g2 <- gene2exprs(rtca,gene2)
##  p <- vector("numeric", length=nrow(g1))
##  p[1] <- p[2] <- 1
##  for(i in 3:nrow(g1)) {
##    if(all(g1[i,]==1) && all(g2[i,]==1)) p[i] <- 1
##    else 
##      p[i] <- t.test(g1[i,], g2[i,], ...)$p.value
##  }
##
##  return(p)
##}

##genesymbolValue <- function(x, gs) {
##  pdata <- pData(x)
##  frame <- exprs(x)
##  gswell <- grep(gs,pdata$GeneSymbol)
##  gscols <- frame[,gswell, drop=FALSE]
##  return(gscols)
##}

## normalization strategy: crop control
##cropCtr <- function(x, field="GeneSymbol", ctrs=c("Allstar")) {
##  pd <- pData(x)
##  stopifnot(field %in% colnames(pd))
##
##  ctrExp <- rowMeans(exprs(x)[,pd[,field] %in% ctrs])
##  exprs(x) <- apply(exprs(x),2, function(y) y-ctrExp)
##  return(x)
##}
##
