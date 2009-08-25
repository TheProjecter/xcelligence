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
