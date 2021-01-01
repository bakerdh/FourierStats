#' tsqc.test: function to calculate the T-squared circ statistic from Victor & Mast (1991)
#' the inputs must be Nx2 or 2xN matrices of numbers
#' for a one-sample test, mu is an optional vector to which the data are compared
#' for a two-sample test, x and y are the data from the two conditions
#' if paired=TRUE, the matrices (x and y) must be the same size
#' @export
tsqc.test <- function(x,y=NULL,paired=FALSE,mu=NULL){

  # transpose if dimensions are the wrong way round
  if (!is.null(ncol(x))){d <- dim(x)
  if (d[1]<d[2]){x <- t(x)}}
  if (!is.null(y)){if (!is.null(ncol(y))){d <- dim(y)
  if (d[1]<d[2]){y <- t(y)}}}

  if (!is.null(mu)){for (n in 1:length(mu)){x[,n] <- x[,n] - mu[n]}}

  # convert data to complex representation
  if (!is.complex(x)){x <- complex(real=x[,1],imaginary=x[,2])}
  if (!is.null(y)){if (!is.complex(y)){y <- complex(real=y[,1],imaginary=y[,2])}}

  if (is.null(y)){method <- 'One-sample t-squared circ test'}
  if (!is.null(y) & paired==TRUE){method <- 'Paired t-squared circ test'}
  if (!is.null(y) & paired==FALSE){method <- 'Independent samples t-squared circ test'}

  if (paired==TRUE){if (is.null(y)==FALSE){x <- x - y}}

  # paired or one-sample version of the test
  if (paired==TRUE | is.null(y)){

    nobs <- length(x)
    cohmean <- mean(x)
    absmean <- abs(cohmean)

    displacement <- abs(x - cohmean)
    diffsumsq <- sum(displacement^2)
    tsqc <- (nobs-1) * (absmean^2) / diffsumsq
    Fratio <- nobs*tsqc
    df2 <- nobs*2-2
    pval <- pf(Fratio,df1=2,df2=df2,lower.tail=FALSE)
  }

  # independent samples version of the test
  if (paired==FALSE & !is.null(y)){
    nobs1 <- length(x)
    cohmean1 <- mean(x)
    absmean1 <- abs(cohmean1)
    displacement1 <- abs(x - cohmean1)

    nobs2 <- length(y)
    cohmean2 <- mean(y)
    absmean2 <- abs(cohmean2)
    displacement2 <- abs(y - cohmean2)

    meandiff <- cohmean1 - cohmean2
    absdiff <- abs(meandiff)
    diffsumsq <- sum(displacement1^2) + sum(displacement2^2)

    tsqc <- (nobs1 + nobs2 - 2) * (absdiff^2) / diffsumsq
    Fratio <- ((nobs1*nobs2)/(nobs1 + nobs2))*tsqc
    df2 <- (2*nobs1 + 2*nobs2 - 4)
    pval <- pf(Fratio,df1=2,df2=df2,lower.tail=FALSE)
  }

  df1 <- 2
  p.value <- min(pval,1)

  output <- data.frame(tsqc,Fratio,df1,df2,p.value,method)
  return(output)}
