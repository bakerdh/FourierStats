#' CI.test: function to calculate the condition index of bivariate data
#' the condition index is the ratio of eigenvalues (eigenvector lengths), calculated longest/shortest
#' the index is then compared for significance with an expected distribution function
#' significant tests (p<0.05) violate the assumptions of the T-squared-circ and ANOVA-squared-circ tests
#' see Baker (2021) for further details
#' @export
CI.test <- function(data){

  if (is.complex(data)){data <- data.frame(Re(data),Im(data))}

  C <- cov(data)  # covariance matrix
  N <- nrow(data) # sample size
  eigVal  <- eigen(C)$values  # calculate eigenvalues
  CI <- sqrt(eigVal[1]/eigVal[2]) # calculate condition index

  maxci <- max(20,round(CI+5))   # make sure our target CI is in range
  cilist <- seq(1,maxci,0.01)    # list of possible CIs

  # calculate probability density function using modified Edelman equation
  pdffunction <- ((N-2)*(2^(N-2))) *
    ((cilist^2 - 1)/((cilist^2+1)^(N-1))) * (cilist^(N-3))
  cdfinverse <- 1-(cumsum(pdffunction)/sum(pdffunction)) # inverse of cdf
  criticalCI <- cilist[min(which(cdfinverse<=0.05))]  # find the threshold CI

  indices <- which(cilist>=CI)  # find the CI values larger than our CI
  pval <- cdfinverse[indices[1]] # estimate the p-value

  output <- data.frame(CI,N,criticalCI,pval)
  return(output)}
