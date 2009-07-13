library(RTCA)


tl <- new("RTCAtimeline")
timeUnit(tl) <- "hour"
startTime(tl) <- Sys.time()-3e6
show(tl)

tl2 <- addAction(tl, 1, "seeding")
show(tl2)
stopifnot(getAction(tl2, 0)=="start")
stopifnot(is.null(getAction(tl2, -1)))

rmAction(tl2, 0)
rmAction(tl2, 1)
rmAction(tl2, -1)

updateAction(tl2, 1, "not seeding")
updateAction(tl2, -1,"not seeding")
updateAction(tl2, -1,"not seeding", ifnotfound="ignore")
startTime(tl2) <- as.POSIXct("2009-05-06 14:52:03")
tl2

ofile <- system.file("/extdata/testOutput.csv", package="RTCA")
x <- parseRTCA(ofile)

## add actions to timeline
x <- addAction(x, 22, "transfection")
x <- addAction(x, 30, "change medium")

xr <- ratioTransform(x, 35)
xrs <- smoothTransform(xr)
xi <- interpolationTransform(x)
xd <- derivativeTransform(x)
xrgr <- rgrTransform(x)
#plotRTCA(xrgr, xlim=c(35, 80), ylim=c(-0.2, 0.2))
plot(sliceRTCA(x,0, 80)[,c(11,13)], type="l", col="black", ylim=c(0, 1.5))
abline(h=0, col=2, lty=5)
plot(sliceRTCA(xrgr, 20, 80)[,c(11,13)], type="l", col="black")
abline(h=0, col=2, lty=5)

plotGridEffect(x)
plotGridEffect(x, "col")
plateView(sliceRTCA(x, 0, 80))
