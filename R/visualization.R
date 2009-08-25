## myheatmap.2: modified the 'heatmap.2' function in the 'gplots' package (2.7.1)
## use 'help(package="gplots")' and 'help("heatmap.2", package=gplots)' to view help
myheatmap.2 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, col = "heat.colors", colsep, 
    rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, 
    notecex = 1, notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE), densadj = 0.25, main = NULL, xlab = NULL, 
    ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if ((Colv == "Rowv") && (!isTRUE(Rowv) || is.null(Rowv))) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) 
        if (missing(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    if (length(breaks) == 1) {
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
            length = breaks)
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    else if (is.character(col) && length(col) == 1) 
        col <- do.call(col, list(ncol))
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[] <- ifelse(x < min.breaks, min.breaks, x)
    x[] <- ifelse(x > max.breaks, max.breaks, x)
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
        x <- t(x)
        cellnote <- t(cellnote)
    }
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      trLine <- apply(x, 2, function(y) peakCenter(rtcawatershed(y)))
      yv <- 1:ncol(x)
      rowFactor <- factor(RowSideColors)
      rowfac <- rev(as.integer(rowFactor))
      for(i in unique(rowfac)) {
        lines(x=trLine[rowfac==i], y=yv[rowfac==i], lwd=2, col=levels(rowFactor)[i], type="s")
    }

    }
    if (trace %in% c("both", "row")) {
        for (i in rowInd) {
            if (!is.null(hline)) {
                hline.vals <- scale01(hline, min.scale, max.scale)
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(4, 2, 2, 1), mgp=c(1,0.5,0),cex = 0.75)
        if (symkey) {
            max.raw <- max(abs(x), na.rm = TRUE)
            min.raw <- -max.raw
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv, tck=-0.1, cex.axis=0.85)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2, cex=0.8)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
##            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
##            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key", cex.main=1)
            par(cex = 0.5)
        }
        else title("Color Key")
    }
    else plot.new()
    invisible(list(rowInd = rowInd, colInd = colInd))
}

##----------------------------------------##
## Plate and gene oriented visualizations
##----------------------------------------##
## special function to determine grid effect
.gridEffect <- function(rtca, mode=c("row","column")) {
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

plotGridEffect <- function(rtca, mode=c("column","row"), xlab="time point", ylab="readout",legend=TRUE,...) {
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
  if(legend) {
    legend("topright",
           legend=1:ncol(mean), lty=1, lwd=4, col=cols, ncol=4, bty="n", title=mode)
  }
}

controlView <- function(rtca,
                        genesymbol=c("Allstar","COPB2","GFP","mock", "PLK1","WEE1"),
                        cols,
                        ylim,smooth=FALSE, group=TRUE, ylab="Normalized cell index", xlab="Time interval (hour)",  drawsd=TRUE, normline=TRUE, ncol=1, legendpos="topleft",...) {
  timeint <- timepoints(rtca)
  frame <- exprs(rtca)
  pdata <- pData(rtca)
  if(missing(ylim))
    ylim <- quantile(frame, c(0.05, 0.95))
  plot(timeint, frame[,1], type="n", ylim=ylim, xlab=xlab, ylab=ylab,...)
  genesymbols <- relevels(factor(genesymbol), genesymbol)
  if(missing(cols))
    cols <- brewer.pal(nlevels(genesymbols),"Set2")
  for (i in 1:nlevels(genesymbols)) {
    gs <- levels(genesymbols)[i]
    gsindex <- grep(paste("^",gs,"$",sep=""),pdata$GeneSymbol)
    if(group) {
      gscols <- frame[,gsindex, drop=FALSE]
      gscols <- gscols[,!apply(gscols, 2, function(x) all(is.na(x))), drop=FALSE]
      groupmean <- rowMeans(gscols)
      groupsd <- apply(gscols,1,sd,na.rm=T)/sqrt(length(gsindex))*1.96
      lines(timeint, groupmean, col=cols[i], lwd=3)
      if(drawsd) {
        segments(timeint, groupmean+groupsd, timeint, groupmean-groupsd, col=cols[i], lwd=1)
        segments(timeint-0.25, groupmean-groupsd, timeint+0.25, groupmean-groupsd, col=cols[i])
        segments(timeint-0.25, groupmean+groupsd, timeint+0.25, groupmean+groupsd, col=cols[i])
      }
    } else {
      for (j in gsindex) {
        y <- frame[,j]
        if(smooth) y <- smooth(y)
        lines(timeint, y, col=cols[i], lwd=3)
      }
    }
  }
  grid()
  legend(legendpos, legend=levels(genesymbols), col=cols, lwd=4, cex=1.1, bty="n", ncol=ncol)


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
  titles <- rownames(pData(rtca))
  for (i in 1:96) {
    plot(timepoints(rtca),express[,i], type="l", xlab="", ylab="", axes=FALSE,ylim=ylim,...)
    title(titles[i])
    abline(h=0, lty=2, col="darkgrey")
  }
  par(opar)
}
