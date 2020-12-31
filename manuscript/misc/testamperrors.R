tic()
N <- 10
simdata <- matrix(rnorm(2*N,mean=0.5,sd=0.2),nrow=N,ncol=2)
simdata <- complex(real=simdata[,1],imaginary=simdata[,2])
qval <- 95
bootstrapped <- amperrors(simdata, method='boot',quantiles=qval)
circle <- amperrors(simdata, method='circ',quantiles=qval)
ellipse <- amperrors(simdata, method='ellipse',quantiles=qval)
absamp <- amperrors(simdata, method='abs',quantiles=qval)
toc()

plot(x=NULL,y=NULL,axes=FALSE,ann=FALSE, xlim=c(0,5), ylim=c(0,1))
ticklocsx <- 0:5   # locations of tick marks on x axis
ticklocsy <- (0:5)/5    # locations of tick marks on y axis
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)

trueCI <- 1.96*(0.2/sqrt(N))
polygon(c(0,5,5,0), 0.5*sqrt(2)+c(-1,-1,1,1)*trueCI, col=rgb(0,0,0,alpha=0.2),border=NA)
lines(c(0,5),0.5*sqrt(c(2,2)),lty=2)


arrows(1,absamp$meanamp,x1=1, y1=absamp$upperCI, length=0.015, angle=90, lwd=2, col='black')
arrows(1,absamp$meanamp,x1=1, y1=absamp$lowerCI, length=0.015, angle=90, lwd=2, col='black')
arrows(2,circle$meanamp,x1=2, y1=circle$upperCI, length=0.015, angle=90, lwd=2, col='black')
arrows(2,circle$meanamp,x1=2, y1=circle$lowerCI, length=0.015, angle=90, lwd=2, col='black')
arrows(3,ellipse$meanamp,x1=3, y1=ellipse$upperCI, length=0.015, angle=90, lwd=2, col='black')
arrows(3,ellipse$meanamp,x1=3, y1=ellipse$lowerCI, length=0.015, angle=90, lwd=2, col='black')
arrows(4,bootstrapped$meanamp,x1=4, y1=bootstrapped$upperCI, length=0.015, angle=90, lwd=2, col='black')
arrows(4,bootstrapped$meanamp,x1=4, y1=bootstrapped$lowerCI, length=0.015, angle=90, lwd=2, col='black')

points(1,absamp$meanamp,pch=21,cex=2,bg='black')
points(2,circle$meanamp,pch=21,cex=2,bg='red')
points(3,ellipse$meanamp,pch=21,cex=2,bg='blue')
points(4,bootstrapped$meanamp,pch=21,cex=2,bg='darkgreen')


legend(0, 0.45, c("Absolute","Circular","Ellipse","Bootstrapped"), cex=0.7, col="black", pt.cex=1.4, pt.bg=c("black","red","blue","darkgreen"), pch=21, pt.lwd=2, box.lwd=2)


