# library(R.matlab)
# library(reshape2)
# data <- readMat('Hwangdata.mat')
# mousedata <- data$mousedata
# montage <- read.csv('~/Desktop/10.12751_g-node.e5tyek/montage.csv')
# outline <- read.csv('~/Desktop/10.12751_g-node.e5tyek/outline.csv',header=FALSE)
# outline <- t(outline[,1:16])
# outline2 <- melt(outline)
# outline <- matrix(outline2[1:268,3],nrow=134,ncol=2,byrow=TRUE)
# save(file='Hwangdata.RData',list=c('mousedata','montage','outline'))

library(HDMD)

v4Interp <- function(df, xo, yo, rmax = .75, gridRes = 67) {
  ## Create a function to perform Matlab's v4 interpolation.
  ## Takes as input a data-frame with columns x, y, and z (x co-ordinates, y co-ordinates, and amplitude)
  ## and variables xo and yo, the co-ordinates which will be use to create a grid for interpolation
  xo <- matrix(rep(xo,length(yo)),nrow = length(xo),ncol = length(yo))
  yo <- t(matrix(rep(yo,length(xo)),nrow = length(yo),ncol = length(xo)))
  xy <- df$x + df$y*sqrt(as.complex(-1))
  d <- matrix(rep(xy,length(xy)),nrow = length(xy), ncol = length(xy))
  d <- abs(d - t(d))
  diag(d) <- 1
  g <- (d^2) * (log(d)-1)   # Green's function.
  diag(g) <- 0
  weights <- qr.solve(g,df$z)
  xy <- t(xy)
  outmat <- matrix(nrow = gridRes,ncol = gridRes)
  for (i in 1:gridRes){
    for (j in 1:gridRes) {
      test4 <- abs((xo[i,j] + sqrt(as.complex(-1))*yo[i,j]) - xy)
      g <- (test4^2) * (log(test4)-1)
      outmat[i,j] <- g %*% weights}}
  outDf <- data.frame(x = xo[,1],outmat)
  names(outDf)[1:length(yo[1,])+1] <- yo[1,]
  return(outDf)}

nbdtpal <- c(rgb(0,0,0),rgb(187/255,40/255,107/255),rgb(56/255,111/255,164/255),'darkgreen','darkorange')

load('Hwangdata.RData')

frequencies <- 1:100
f1 <- which(frequencies==40)

meandata <- abs(apply(mousedata,c(2,3,4),mean))

temp <- rowMeans(mousedata[,1:2,11,f1])
xdata <- Re(temp)
ydata <- Im(temp)
temp2 <- rowMeans(mousedata[,1:2,12,f1])
xdata2 <- Re(temp2)
ydata2 <- Im(temp2)


plot(x=NULL,y=NULL,axes=FALSE,ann=FALSE, xlim=c(0,150), ylim=c(0,2))
ticklocsx <- seq(0,80,20)    # locations of tick marks on x axis
ticklocsy <- seq(0,2,0.5)    # locations of tick marks on y axis
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklocsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklocsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
title(ylab="Amplitude (µV)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
mtext('Frequency (Hz)',side=1,at=40,line=1.2,cex=1.5)

rmax <- 7   #specify a maximum boundary for the grid
gridRes <- 100 #specify the interpolation grid resolution
ramp <- colorRamp(c("white", rgb(187/255,40/255,107/255)))  # create a ramp from one colour to another
colmatrix2 <- rgb(ramp(seq(0, 1, length = 100)), max = 255)   # index the ramp at ten points

datatoplot <- abs(meandata[,11,f1])
testDat<- data.frame(x = montage[,2],
                     y = -montage[,3],
                     z = datatoplot)

xo <- seq(min(-rmax, testDat$x), max(rmax, testDat$x), length = gridRes)
yo <- seq(max(rmax, testDat$y), min(-rmax, testDat$y), length = gridRes)

interpV4 <- v4Interp(testDat, xo, yo, rmax, gridRes)

zo2 <- as.matrix(interpV4[,2:ncol(interpV4)])
xo2 <- matrix(rep(xo,length(yo)),nrow = length(xo),ncol = length(yo))
yo2 <- t(matrix(rep(yo,length(xo)),nrow = length(yo),ncol = length(xo)))
zo2[which(zo2>2)] <- 2

image(2.4*xo+20,1.5+xo/14,zo2,zlim=c(0,2),col=colmatrix2,add=TRUE,useRaster=TRUE)

rmax2 <- rmax+0.1
maskx <- outline[which(outline[,1]>=0),1]
masky <- outline[which(outline[,1]>=0),2]
polygon(2.4*c(0,maskx,0,rmax2,rmax2,0)+20,1.5+c(rmax2,masky,-rmax2,-rmax2,rmax2,rmax2)/14,border=NA,col="white")
maskx <- outline[which(outline[,1]<=0),1]
masky <- outline[which(outline[,1]<=0),2]
i <- sort(masky,index.return=TRUE)
maskx <- maskx[i$ix]
masky <- masky[i$ix]
polygon(2.4*c(0,maskx,0,-rmax2,-rmax2,0)+20,1.5+c(-rmax2,masky,rmax2,rmax2,-rmax2,-rmax2)/14,border=NA,col="white")

lines(2.4*outline[,1]+20,1.5+outline[,2]/14,lwd=2)
points(2.4*montage[3:38,2]+20,1.5+montage[3:38,3]/14,pch=16,col='grey')
points(2.4*montage[1:2,2]+20,1.5+montage[1:2,3]/14,pch=16)


ramp <- colorRamp(c("white", rgb(56/255,111/255,164/255)))  # create a ramp from one colour to another
colmatrix2 <- rgb(ramp(seq(0, 1, length = 100)), max = 255)   # index the ramp at ten points

datatoplot <- abs(meandata[,12,f1])
testDat<- data.frame(x = montage[,2],
                     y = -montage[,3],
                     z = datatoplot)

xo <- seq(min(-rmax, testDat$x), max(rmax, testDat$x), length = gridRes)
yo <- seq(max(rmax, testDat$y), min(-rmax, testDat$y), length = gridRes)

interpV4 <- v4Interp(testDat, xo, yo, rmax, gridRes)

zo2 <- as.matrix(interpV4[,2:ncol(interpV4)])
xo2 <- matrix(rep(xo,length(yo)),nrow = length(xo),ncol = length(yo))
yo2 <- t(matrix(rep(yo,length(xo)),nrow = length(yo),ncol = length(xo)))
zo2[which(zo2>2)] <- 2

image(2.4*xo+60,1.5+xo/14,zo2,zlim=c(0,2),col=colmatrix2,add=TRUE,useRaster=TRUE)

rmax2 <- rmax+0.1
maskx <- outline[which(outline[,1]>=0),1]
masky <- outline[which(outline[,1]>=0),2]
polygon(2.4*c(0,maskx,0,rmax2,rmax2,0)+60,1.5+c(rmax2,masky,-rmax2,-rmax2,rmax2,rmax2)/14,border=NA,col="white")
maskx <- outline[which(outline[,1]<=0),1]
masky <- outline[which(outline[,1]<=0),2]
i <- sort(masky,index.return=TRUE)
maskx <- maskx[i$ix]
masky <- masky[i$ix]
polygon(2.4*c(0,maskx,0,-rmax2,-rmax2,0)+60,1.5+c(-rmax2,masky,rmax2,rmax2,-rmax2,-rmax2)/14,border=NA,col="white")

lines(2.4*outline[,1]+60,1.5+outline[,2]/14,lwd=2)
points(2.4*montage[3:38,2]+60,1.5+montage[3:38,3]/14,pch=16,col='grey')
points(2.4*montage[1:2,2]+60,1.5+montage[1:2,3]/14,pch=16)



lines(frequencies[1:79],meandata[1,11,1:79],col=nbdtpal[2],lwd=3)
lines(frequencies[1:79],meandata[1,12,1:79],col=nbdtpal[3],lwd=3)
lines(c(90,150),c(1,1),lwd=2)
lines(c(120,120),c(0,2),lwd=2)

ellRot1 <- getel(data.frame(xdata,ydata))
ellRot2 <- getel(data.frame(xdata2,ydata2))
polygon(10*(ellRot1[1, ]+mean(xdata))+120, 1+(ellRot1[2, ]+mean(ydata))/3, col=addalpha(nbdtpal[2],0.3),border=NA)
polygon(10*(ellRot2[1, ]+mean(xdata2))+120, 1+(ellRot2[2, ]+mean(ydata2))/3, col=addalpha(nbdtpal[3],0.3),border=NA)

points(10*xdata+120,1+ydata/3,pch=16,col=nbdtpal[2])
points(10*xdata2+120,1+ydata2/3,pch=16,col=nbdtpal[3])
points(10*mean(xdata)+120,1+mean(ydata)/3,pch=21,bg=nbdtpal[2],cex=1.5)
points(10*mean(xdata2)+120,1+mean(ydata2)/3,pch=21,bg=nbdtpal[3],cex=1.5)

legend(125,0.8,c('Sound','Light'),pch=21,pt.cex=1.5,pt.bg=nbdtpal[2:3],box.lwd=2)

text(20,1.92,'Sound',cex=1.5,adj=0.5)
text(60,1.92,'Light',cex=1.5,adj=0.5)

text(0,1.95,'(a)',cex=2,adj=0.5)
text(90,1.95,'(b)',cex=2,adj=0.5)

text(149,0.95,'3',cex=1.5,adj=0.5)
text(121.5,1.97,'3',cex=1.5,adj=0.5)



totest <- data.frame(xdata,ydata)
totest2 <- data.frame(xdata2,ydata2)
CI.test(totest)
CI.test(totest2)
tsqc.test(totest,y=totest2,paired=TRUE)

abs1 <- abs(rowMeans(mousedata[,1:2,11,f1]))
abs2 <- abs(rowMeans(mousedata[,1:2,12,f1]))
t.test(abs1,abs2,paired=TRUE)

colnames(totest2) <- colnames(totest)
fores <- rbind(totest,totest2)
mahala_sq <- pairwise.mahalanobis(fores, rep(1:2,each=6))
D <- sqrt(mahala_sq$distance[1,2])

