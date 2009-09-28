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



