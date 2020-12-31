#' amperrors: calculate error bars on the amplitude component of coherently averaged data
#' the input is a vector of complex numbers, or an Nx2 matrix of real and imaginary components
#' the optional 'method' flag allows the following values:
#'  - 'boot': (default) bootstraps the confidence intervals (10000 iterations)
#'  - 'circ': calculates confidence intervals based on a circular bounding region (used when the condition index test is non-significant)
#'  - 'ellipse': calculates confidence intervals based on an elliptical bounding region (based on Pei et al. 2017, doi: 10.1016/j.visres.2016.03.010)
#'  - 'abs': calculates the standard error (and the mean) using the absolute amplitudes
#'  the optional 'quantiles' flag allows the quantile to be set for the confidence intervals
#'  typical values are 95 (for 95% confidence intervals), and 68 (for standard errors)
#'  the function returns the mean amplitude and the upper and lower error bars
#' @export
amperrors <- function(input, method='boot',quantiles=95){
  d <- dim(input)
  compdata <- input
  if (!is.complex(input)){compdata <- complex(real=input[,1],imaginary=input[,2])}
  xydata <- data.frame(Re(compdata),Im(compdata))
  # compdata is now expressed in complex values, xy data is an Nx2 data frame of real and imaginary values
  output <- NULL
  output$meanamp <- abs(mean(compdata)) # calculate the mean amplitude

  if (method=='boot'){
    # bootstrap the error bars on the amplitude using 10000 resamples
    bspop <- NULL
    for (n in 1:10000){bspop[n] <- abs(mean(sample(compdata,length(compdata),replace=TRUE)))}
    output$upperCI <- quantile(bspop,probs=1-(1-(quantiles/100))/2)
    output$lowerCI <- quantile(bspop,probs=(1-(quantiles/100))/2)
  }

  if (method=='circ'){
    # pools the variance in the x (real) and y (imaginary) directions
    # calculates the standard error based on this estimate
    sdcirc <- sqrt((sd(xydata[,2])^2 + sd(xydata[,2])^2)/2)
    se <- sdcirc/sqrt(length(compdata))
    if (quantiles==95){ci <- 1.96*se}
    if (quantiles==68){ci <- se}
    output$upperCI <- output$meanamp + ci
    output$lowerCI <- output$meanamp - ci
  }

  if (method=='ellipse'){
    # implements the method of Pei et al. (2017)
    # calculates the bounding ellipse for the data points
    # then uses this to find the nearest and farthest points from the origin
    # these points are used as the lower and upper error bars

    angles <- seq(0, 2*pi, length.out=200)          # angles for ellipse
    evs <- eigen(cov(xydata))
    eigVal  <- 1.96*eigen(cov(xydata))$values/sqrt(nrow(xydata))
    eigVec  <- eigen(cov(xydata))$vectors
    eigScl  <- eigVec  %*% diag(sqrt(eigVal))
    xMat    <- rbind(Re(mean(compdata)) + eigScl[1, ], Re(mean(compdata)) - eigScl[1, ])
    yMat    <- rbind(Im(mean(compdata)) + eigScl[2, ], Im(mean(compdata)) - eigScl[2, ])
    ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles))
    ellRot  <- eigVec %*% t(ellBase)
    ellRot[1,] <- ellRot[1,]+Re(mean(compdata))
    ellRot[2,] <- ellRot[2,]+Im(mean(compdata))

    totaldistfrom0 <- sqrt(ellRot[1,]^2 + ellRot[2,]^2)
    i <- which(totaldistfrom0==min(totaldistfrom0))
    j <- which(totaldistfrom0==max(totaldistfrom0))
    ellCI <- abs(complex(real=ellRot[1,c(i,j)],imaginary=ellRot[2,c(i,j)]))
    output$upperCI <- ellCI[1]
    output$lowerCI <- ellCI[2]
    if (quantiles==68){output$upperCI <- ((ellCI[1] - output$meanamp)/1.96) + output$meanamp}
    if (quantiles==68){output$lowerCI <- ((ellCI[2] - output$meanamp)/1.96) + output$meanamp}
  }

  if (method=='abs'){
    # calculate the mean and the error bars using absolute amplitude values (incoherent averaging)
    adata <- abs(compdata)
    output$meanamp <- mean(adata) # calculate the mean absolute amplitude
    se <- sd(adata)/sqrt(length(adata))
    if (quantiles==95){ci <- 1.96*se}
    if (quantiles==68){ci <- se}
    output$upperCI <- output$meanamp + ci
    output$lowerCI <- output$meanamp - ci
  }

return(output)}
