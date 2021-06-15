
# install the FourierStats package from github using the devtools package
packagelist <- c('FourierStats')
missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
if (length(missingpackages)>0){devtools::install_github("bakerdh/FourierStats")}
toinstall <- packagelist[which(!packagelist %in% (.packages()))]
invisible(lapply(toinstall,library,character.only=TRUE))
addalpha <- function(col, alpha=1){apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))}

library(tictoc)

tic()

nsims <- 1000000

samplesizes <- 2^(2:6)
nconds <- c(3,5,7,9)
finex <- seq(0,4,length=1000)


pdf('Fdistributions.pdf', bg="transparent", height = 5, width = 9)

plotlims <- c(0,30,0,4)
ticklocsx <- seq(0,4,1)    # locations of tick marks on x axis
ticklocsy <- 0:4    # locations of tick marks on y axis
plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# mtext(text = ticklocsy[1:2], side = 2, at=ticklocsy[1:2], line=0.2, las=1)
title(xlab="F ratio       ", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
title(ylab="Probability density", line=1, cex.lab=1.5)

Fmatrix <- array(0,dim=c(length(samplesizes),length(nconds),nsims))

for (s in 1:length(samplesizes)){

axis(1, at=ticklocsx+(6*(s-1)), tck=0.01, lab=F, lwd=2)
mtext(text = ticklocsx, side = 1, line=0.2, at=ticklocsx+(6*(s-1)))
mtext(text = paste('N =',samplesizes[s]), side=3, at=2+(6*(s-1)))

       for (s2 in 1:length(nconds)){

      grouplabels <- as.factor(rep(1:nconds[s2],each=samplesizes[s]))

      for (n in 1:nsims){

      data <- matrix(rnorm(nconds[s2]*2*samplesizes[s]),nrow=nconds[s2]*samplesizes[s],ncol=2)

      simdata <- complex(real=data[,1],imaginary=data[,2])

      output <- anovacirc.test(data.frame(simdata,grouplabels))
      Fmatrix[s,s2,n] <- output$Fratio

}



a <- hist(Fmatrix[s,s2,],breaks = 200, plot=FALSE)
axvals <- a$mids
ayvals <- 0.95*a$counts/max(a$counts)
polygon(c(0,axvals,max(axvals))+(6*(s-1)), c(0,ayvals,0)+(4-s2), col=addalpha('blue',alpha=0.2),border=NA)
dfM <- 2*(nconds[s2]-1)
dfR <- 2*((samplesizes[s]*nconds[s2]) - nconds[s2])
fdist <- df(finex,dfM,dfR)
fdist <- 0.95*fdist/max(fdist)
lines(finex+(6*(s-1)),(4-s2)+fdist,lwd=3)

}}

text(29,3.5,'k = 3')
text(29,2.5,'k = 5')
text(29,1.5,'k = 7')
text(29,0.5,'k = 9')

dev.off()

toc()

