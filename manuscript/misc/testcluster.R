
N <- 10
m <- 12

data <- matrix(rnorm(N*m),nrow=N,ncol=12)
data[,3:6] <- data[,3:6] + 1

allp <- NULL
for (i in 1:m){
  output <- t.test(data[,i])
  allp[i] <- output$p.value
}

clusts <- clustercorrect(data,testtype=1,clustformthresh=0.05,nresamples=1000)



load("/Users/danbaker/Google Drive/Current work/Research/powercontours/headplotspec.RData")

load('~/Google Drive/Current work/R scripts/danlabtoolbox/data/montagedata.RData')

distmatrix <- matrix(0,nrow=64,ncol=64)
adjacencymatrix <- matrix(0,nrow=64,ncol=64)
for (n in 1:64){
  for (m in 1:64){
    xy1 <- montage$electrodelocs[montage$channelmappings[n],1:2]
    xy2 <- montage$electrodelocs[montage$channelmappings[m],1:2]
    distmatrix[n,m] <- as.numeric(sqrt((xy1[1]-xy2[1])^2 + (xy1[2]-xy2[2])^2))
  }
}
threshold <- 0.18
i1 <- which(distmatrix<=threshold)
i2 <- which(distmatrix==0)
adjacencymatrix[i1] <- 1
adjacencymatrix[i2] <- 0

par(pty="s")  # make axis square
plot(x=NULL,y=NULL,axes=FALSE,ann=FALSE, xlim=c(-0.6,0.6), ylim=c(-0.6,0.6))

for (n in 1:64){
  for (m in 1:64){
    if (adjacencymatrix[n,m]){
      lines(montage$electrodelocs[montage$channelmappings[c(n,m)],1],montage$electrodelocs[montage$channelmappings[c(n,m)],2])
    }
}}
points(montage$electrodelocs[1:64,1],montage$electrodelocs[1:64,2],pch=16,col='grey')

lines(montage$headoutline[,1],montage$headoutline[,2],col="black",lwd=2)
lines(montage$noseoutline[,1],montage$noseoutline[,2],col="black",lwd=2)
lines(montage$Rearoutline[,1],montage$Rearoutline[,2],col="black",lwd=2)
lines(montage$Learoutline[,1],montage$Learoutline[,2],col="black",lwd=2)





ramp <- colorRamp(c("white", rgb(56/255,111/255,164/255)))
colmatrix2 <- rgb(ramp(seq(0, 1, length = 100)), max = 255)   # index the 

cond <- 2
datatoplot <- abs(apply(alltarget[,1:64,cond,], 2, mean))
datatoplot[c(13,19)] <- mean(datatoplot)   # set the mastoids to 0
testDat<- data.frame(x = montage$electrodelocs[montage$channelmappings,1],
                     y = -montage$electrodelocs[montage$channelmappings,2],
                     z = datatoplot)

rmax <- 0.55   #specify a maximum boundary for the grid
gridRes <- 100 #specify the interpolation grid resolution

xo <- seq(min(-rmax, testDat$x), max(rmax, testDat$x), length = gridRes)
yo <- seq(max(rmax, testDat$y), min(-rmax, testDat$y), length = gridRes)

interpV4 <- v4Interp(testDat, xo, yo, rmax, gridRes)

zo2 <- as.matrix(interpV4[,2:ncol(interpV4)])
xo2 <- matrix(rep(xo,length(yo)),nrow = length(xo),ncol = length(yo))
yo2 <- t(matrix(rep(yo,length(xo)),nrow = length(yo),ncol = length(xo)))
outsidecircle <- sqrt(xo2^2 + yo2^2) > 0.5
zo2[outsidecircle] <- 0

plot(x=NULL,y=NULL,axes=FALSE,ann=FALSE, xlim=c(-0.6,0.6), ylim=c(-0.6,0.6))

image(xo,xo,zo2,zlim=c(0,0.45),col=colmatrix2,add=TRUE,useRaster=TRUE)

for (n in 1:64){points(montage$electrodelocs[n,1],montage$electrodelocs[n,2],pch=16,col="grey")}

lines(montage$headoutline[,1],montage$headoutline[,2],col="black",lwd=2)
lines(montage$noseoutline[,1],montage$noseoutline[,2],col="black",lwd=2)
lines(montage$Rearoutline[,1],montage$Rearoutline[,2],col="black",lwd=2)
lines(montage$Learoutline[,1],montage$Learoutline[,2],col="black",lwd=2)


dataforstats <- apply(alltarget[,1:64,cond,],1:2,mean)

sigpoints <- NULL
for (n in 1:64){
  temp <- dataforstats[,n]
  output <- tsqc.test(temp)
  if (output$p.value<0.05){sigpoints <- c(sigpoints,n)}
}

# meandata <- abs(colMeans(dataforstats))
# i <- which(meandata==max(meandata))
points(montage$electrodelocs[montage$channelmappings[sigpoints],1],montage$electrodelocs[montage$channelmappings[sigpoints],2],pch=16,col="red")

clusts <- clustercorrect(dataforstats,testtype=3,adjacencymatrix=adjacencymatrix,clustformthresh=0.05,nresamples=1000)

toplot <- clusts$clusterpoints[[1]]
points(montage$electrodelocs[montage$channelmappings[toplot],1],montage$electrodelocs[montage$channelmappings[toplot],2],pch=16,col="red")













load('FourierStats/manuscript/misc/Hwangdata.RData')

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

distmatrix <- matrix(0,nrow=38,ncol=38)
adjacencymatrix <- matrix(0,nrow=38,ncol=38)
for (n in 1:38){
  for (m in 1:38){
    xy1 <- montage[n,2:3]
    xy2 <- montage[m,2:3]
    distmatrix[n,m] <- as.numeric(sqrt((xy1[1]-xy2[1])^2 + (xy1[2]-xy2[2])^2))
  }
}
threshold <- 1.8
i1 <- which(distmatrix<=threshold)
i2 <- which(distmatrix==0)
adjacencymatrix[i1] <- 1
adjacencymatrix[i2] <- 0

# par(pty="s")  # make axis square
# plot(x=NULL,y=NULL,axes=FALSE,ann=FALSE, xlim=c(-7,7), ylim=c(-7,7))
# lines(outline[,1],outline[,2],lwd=2)
# for (n in 1:38){
#   for (m in 1:38){
#     if (adjacencymatrix[n,m]){
#       lines(montage[c(n,m),2],montage[c(n,m),3])
#     }
# }}
# points(montage[1:38,2],montage[1:38,3],pch=16,col='grey')




frequencies <- 1:100
f1 <- which(frequencies==40)

meandata <- abs(apply(mousedata,c(2,3,4),mean))

temp <- rowMeans(mousedata[,1:2,11,f1])
xdata <- Re(temp)
ydata <- Im(temp)
temp2 <- rowMeans(mousedata[,1:2,12,f1])
xdata2 <- Re(temp2)
ydata2 <- Im(temp2)

thiscond <- 12

par(pty="s")  # make axis square
plot(x=NULL,y=NULL,axes=FALSE,ann=FALSE, xlim=c(-7,7), ylim=c(-7,7))

rmax <- 7   #specify a maximum boundary for the grid
gridRes <- 100 #specify the interpolation grid resolution
ramp <- colorRamp(c("white", rgb(0.22,0.44,0.64)))  # create a ramp from one colour to another
colmatrix2 <- rgb(ramp(seq(0, 1, length = 100)), max = 255)   # index the ramp at ten points

datatoplot <- abs(meandata[,thiscond,f1])
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

image(xo,xo,zo2,zlim=c(0,2),col=colmatrix2,add=TRUE,useRaster=TRUE)

rmax2 <- rmax+0.1
maskx <- outline[which(outline[,1]>=0),1]
masky <- outline[which(outline[,1]>=0),2]
polygon(c(0,maskx,0,rmax2,rmax2,0),c(rmax2,masky,-rmax2,-rmax2,rmax2,rmax2),border=NA,col="white")
maskx <- outline[which(outline[,1]<=0),1]
masky <- outline[which(outline[,1]<=0),2]
i <- sort(masky,index.return=TRUE)
maskx <- maskx[i$ix]
masky <- masky[i$ix]
polygon(c(0,maskx,0,-rmax2,-rmax2,0),c(-rmax2,masky,rmax2,rmax2,-rmax2,-rmax2),border=NA,col="white")

dataforcc <- mousedata[,,thiscond,f1]
clusts <- clustercorrect(dataforcc,adjacencymatrix=adjacencymatrix,testtype=3,nresamples=10000)


lines(outline[,1],outline[,2],lwd=2)
points(montage[,2],montage[,3],pch=16,col='grey')
for (n in 1:clusts$nclusters){
  temp <- clusts$clusterpoints[[n]]
  points(montage[temp,2],montage[temp,3],pch=16,col='red')}

