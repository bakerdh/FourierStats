nbdtpal <- c(rgb(0,0,0),rgb(0.73,0.16,0.42),rgb(0.22,0.44,0.64),'darkgreen','darkorange')
addalpha <- function(col, alpha=1){apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))}
set.seed(170382)

getbanana <- function(p){
  thetaRange <- seq(0.01, pi, 0.01)
  Tval <- atan((p[3]*tan(thetaRange))/p[4])
  x1 <- p[1] + (p[3] * cos(Tval))
  x2 <- p[1] - (p[3] * cos(Tval))
  y <- p[2] + (p[4] * sin(Tval))
  
  XXX <- c(x1[1:(length(x1)/2)], x2[(length(x2)/2):1])
  Y <- c(y[1:(length(y)/2)], y[(length(y)/2):1], y[length(y):(1+length(y)/2)], y[(1+length(y)/2):length(y)])
  X <- c(XXX, XXX[length(XXX):1])
  errormatrix <- pol2cart(matrix(data=c(Y,X),ncol=2))
  return(errormatrix)
}

par(mfrow=c(2,1), mar=c(3,4,1,2))

plotlims <- c(0,2.4,0,6) 
ticklocsx <- seq(0,1,0.5)    # locations of tick marks on x axis
ticklocsy <- seq(0,0.1,0.02)    # locations of tick marks on y axis
ticklabelsx <- ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   
axis(1, at=c(0,1), tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
# axis(2, at=c(0,6), tck=0.01, lab=F, lwd=2)
axis(1, at=c(1.4,2.4), tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
# axis(4, at=c(0,6), tck=0.01, lab=F, lwd=2)
mtext(text = 'Time', side = 1, cex=1.5, line=0.5, at=0.5)     # add the tick labels
mtext(text = 'Frequency', side = 1, cex=1.5, line=0.5, at=1.9)     # add the tick labels

text(1.2,5.5,'Signal',cex=1.2,adj=0.5)
text(1.2,4.5,'Observation',cex=1,adj=0.5)
text(1.2,3.75,'1',cex=1.2,adj=0.5)
text(1.2,2.75,'2',cex=1.2,adj=0.5)
text(1.2,1.75,'3',cex=1.2,adj=0.5)
text(1.2,1.1,'.',cex=1.2,adj=0.5)
text(1.2,1,'.',cex=1.2,adj=0.5)
text(1.2,0.9,'.',cex=1.2,adj=0.5)
text(1.2,0.25,'N',cex=1.2,adj=0.5)

text(0,5.8,'(a)',cex=1.5,adj=0.5)
text(1.45,5.8,'(b)',cex=1.5,adj=0.5)


tvals <- seq(0,10-0.001,0.001)
fvals <- c(0:(length(tvals)/2-1),(length(tvals)/2-1):0)/10
signal <- 0.8*cos(5*tvals * 2*pi + pi)

lines(tvals[1:1000],5+(signal[1:1000]+1)/2,lwd=2)

fsignal <- abs(fft(signal))/length(signal)
lines(1.4+fvals[10:200]/20,5+(2*fsignal[10:200]),lwd=2)

signal <- 0.8*cos(5*tvals * 2*pi + 3*pi/2)
ypos <- c(0,1.5,2.5,3.5)
for (n in 1:4){
nspec <- ((fft(signal)/length(signal)) + 0.2*rnorm(length(fvals))/fvals)
nspec[c(1:10,9991:10000)] <- 0
lines(1.4+fvals[10:200]/20,ypos[n]+(2*abs(nspec[10:200])),lwd=2)

nsignal <- Re(fft(nspec, inverse=TRUE))
nsignal <- nsignal[1:1000]/max(abs(nsignal[1:1000]))
nsignal <- nsignal - mean(nsignal)
lines(tvals[1:1000],ypos[n]-0.5+(nsignal[1:1000]+1),lwd=2)
}

signal <- 0.4*cos(5*tvals * 2*pi + pi/4)
allsimsignals <- NULL
nspec <- fft(signal)/length(signal)
for (n in 1:20){
  respec <- Re(nspec) + 1*rnorm(length(fvals))/fvals
  imspec <- Im(nspec) + 1*rnorm(length(fvals))/fvals
  
  allsimsignals[n] <- complex(real=respec[51],imaginary=imspec[51])
}


plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=c(-1,3.8), ylim=c(-1,1))   

rangles <- seq(0,360)*(pi/180)
r <- seq(0.2,1,0.2)
for (n in 1:length(r)){
  lines(r[n]*cos(rangles),r[n]*sin(rangles),col='gray')
}
lines(c(0,0),c(-1,1),lwd=2)
lines(cos(c(30,210)*pi/180),sin(c(30,210)*pi/180),lwd=2)
lines(cos(c(150,330)*pi/180),sin(c(150,330)*pi/180),lwd=2)

text(0.5,-0.2,'Amplitude',cex=1.2,adj=0.5,srt=-30)
text(1.1,0,'Phase',cex=1.2,adj=0.5,srt=-90)

stdresps <-sd(Mod(allsimsignals))  #/sqrt(length(allsimsignals))
stdangles <- sd(Arg(allsimsignals)) #/sqrt(length(allsimsignals))
pvals <- getbanana(c(mean(Mod(allsimsignals)),Arg(mean(allsimsignals)),stdresps,stdangles))
polygon(pvals[,1], pvals[,2], col=addalpha(nbdtpal[3],0.2), border=NA) # plot upside down hanging from the top axis with our transparent colour

for (n in 1:length(allsimsignals)){lines(c(0,Re(allsimsignals[n])),c(0,Im(allsimsignals[n])),col=addalpha(nbdtpal[1],0.5))}
points(Re(allsimsignals),Im(allsimsignals),pch=16,col=addalpha(nbdtpal[1],0.5),cex=0.8)

meanangle <- mean(((Arg(allsimsignals)-1) %% 2*pi)+1)
meanangle <- Arg(mean(allsimsignals))
meanamp <- mean(Mod(allsimsignals))
y <- meanamp*sin(meanangle)
x <- meanamp*cos(meanangle)
lines(c(0,x),c(0,y),lwd=3,col=nbdtpal[3])
points(x,y,pch=21,bg=nbdtpal[3],cex=1.5)



arrows(2.8,-1,2.8,1,length=0.1,lwd=2)
arrows(1.8,0,3.8,0,length=0.1,lwd=2)
text(1.93,0.1,'Real',cex=1.2,adj=0.5)
text(2.7,-0.7,'Imaginary',cex=1.2,adj=0.5,srt=90)

points(Re(allsimsignals)+2.8,Im(allsimsignals),pch=16,col=addalpha(nbdtpal[1],0.5),cex=0.8)
points(Re(mean(allsimsignals))+2.8,Im(mean(allsimsignals)),pch=21,bg=nbdtpal[3],cex=1.5)
eldata <- getel(allsimsignals)
polygon(eldata[1,]+2.8,eldata[2,], col=addalpha(nbdtpal[3],0.2), border=NA) # plot upside down hanging from the top axis with our transparent colour


text(-1,1,'(c)',cex=1.5,adj=0.5)
text(1.8,1,'(d)',cex=1.5,adj=0.5)

