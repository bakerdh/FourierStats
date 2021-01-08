#' tsqh.test: function to calculate a different variants of the T-squared test of Hotelling (1931)
#' the inputs must be Nx2 or 2xN matrices of numbers, or vectors of complex numbers
#' for a one-sample test, mu is an optional vector to which the data are compared
#' for a two-sample test, x and y are the data from the two conditions
#' if paired=TRUE, the inputs (x and y) must be the same size
#' @export
tsqh.test <- function(x,y=NULL,paired=FALSE,mu=NULL){

  # transpose if dimensions are the wrong way round
  if (!is.null(ncol(x))){d <- dim(x)
  if (d[1]<d[2]){x <- t(x)}}
  if (!is.null(y)){if (!is.null(ncol(y))){d <- dim(y)
  if (d[1]<d[2]){y <- t(y)}}}

  if (!is.null(mu)){for (n in 1:length(mu)){x[,n] <- x[,n] - mu[n]}}

  # convert data to bivariate representation if complex
  if (is.complex(x)){x <- data.frame(Re(x),Im(x))}
  if (!is.null(y)){if (is.complex(y)){y <- data.frame(Re(y),Im(y))}}

  if (is.null(y)){method <- 'One-sample T-squared test'}
  if (!is.null(y) & paired==TRUE){method <- 'Paired T-squared test'}
  if (!is.null(y) & paired==FALSE){method <- 'Independent samples T-squared test'}

  if (paired==TRUE){if (is.null(y)==FALSE){x <- x - y}}

  # paired or one-sample version of the test
  if (paired==TRUE | is.null(y)){

    xbar <- colMeans(x) # estimate means for each DV
    m <- ncol(x)   # record number of DVs
    N <- nrow(x)   # record sample size
    C <- cov(x)    # calculate covariance matrix
    Cinv <- solve(C)  # calculate inverse matrix
    tsq <- N*t(xbar)%*%Cinv%*%xbar # compute t2 statistic
    Fratio <- tsq * (N-m)/(m*(N-1)) # convert to F-ratio
    df1 <- m      # record degrees of freedom
    df2 <- N-m
    pval <- pf(Fratio,df1,df2,lower.tail=FALSE) # estimate p-value

  }

  # independent samples version of the test
  if (paired==FALSE & !is.null(y)){

    xbarA <- colMeans(x) # estimate means for each DV
    m <- ncol(x)   # record number of DVs
    N1 <- nrow(x)   # record sample size

    xbarB <- colMeans(y) # estimate means for each DV
    N2 <- nrow(y)   # record sample size

    CA <- cov(x)    # calculate covariance matrices
    CB <- cov(y)

    meandiff <- xbarA - xbarB

    # combine covariance matrices
    C <- ((N1-1)*CA + (N2-1)*CB)/(N1+N2-m)
    Cinv <- solve(C)                      # calculate inverse matrix

    tsq <- ((N1*N2)/(N1+N2))*t(meandiff)%*%Cinv%*%meandiff    # compute t2 statistic
    Fratio <- tsq * (N1+N2-m-1)/(m*(N1+N2-2))                 # convert to F-ratio
    df1 <- m                                                  # record degrees of freedom
    df2 <- N1+N2-m-1
    pval <- pf(Fratio,df1,df2,lower.tail=FALSE) # estimate p-value

  }

  pval <- min(pval,1)

  # store outputs in a data structure
  output <- data.frame(tsq,Fratio,df1,df2,pval,method)

  return(output)}
