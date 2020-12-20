tsq1.test <- function(data,mu=NULL){
  # function to calculate a one-sample t-squared test after Hotelling (1931)
  # for two-sample tests, use the Hotelling package
  
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
