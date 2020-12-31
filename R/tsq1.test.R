#' tsq1.test: function to calculate a one-sample T-squared test of Hotelling (1931)
#' the input is an N (observations) by m (dependent variables) matrix
#' optional input mu is a vector with m entries, to which the data are compared (otherwise compared to 0)
#' a repeated measures version can be calculated by subtracting the two conditions and then running the one-sample test
#' for two-sample tests, use the hotelling.test function from the Hotelling package
#' @export
tsq1.test <- function(data,mu=NULL){

  # if data are complex, convert to real and absolute values
  if (is.complex(data)){data <- data.frame(Re(data),Im(data))}
  if (is.complex(mu)){mu <- c(Re(mu),Im(mu))}
  # compare each DV to some value
  if (!is.null(mu)){for (m in 1:length(mu)){data[,m] <- data[,m]-mu[m]}}

  xbar <- colMeans(data) # estimate means for each DV
  m <- ncol(data)   # record number of DVs
  N <- nrow(data)   # record sample size
  C <- cov(data)    # calculate covariance matrix
  Cinv <- solve(C)  # calculate inverse matrix
  tsq <- N*t(xbar)%*%Cinv%*%xbar # compute t2 statistic
  Fratio <- tsq * (N-m)/(m*(N-1)) # convert to F-ratio
  df1 <- m      # record degrees of freedom
  df2 <- N-m
  pval <- 1- pf(Fratio,df1,df2) # estimate p-value
  p.value <- min(pval,1)

  # store outputs in a data structure
  output <- data.frame(tsq,Fratio,df1,df2,p.value)
  return(output)}
