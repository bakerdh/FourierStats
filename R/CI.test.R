#' CI.test: function to calculate the condition index of bivariate data
#' Inputs:
#'   data is an N x 2 matrix of bivariate (x,y) observations, or a vector of N complex values
#'   alpha is the criterion for statistical significance, set to a default of 0.05
#' the condition index is square root of the the ratio of eigenvalues (eigenvector lengths), calculated longest/shortest
#' the index is then compared for significance with an expected distribution function
#' significant tests (p < alpha) violate the assumptions of the T-squared-circ and ANOVA-squared-circ tests
#' see Baker (2021) for further details
#' @export
CI.test <- function(data, alpha=0.05){

  if (is.complex(data)){data <- data.frame(Re(data),Im(data))}

  C <- cov(data)  # covariance matrix
  N <- nrow(data) # sample size
  eigVal  <- eigen(C)$values  # calculate eigenvalues
  CI <- sqrt(eigVal[1]/eigVal[2]) # calculate condition index

  cilist <- seq(1,100,0.001)    # list of possible CIs

  # calculate probability density function using modified Edelman equation
  pdffunction <- ((N-2)*(2^(N-2))) *
    ((cilist^2 - 1)/((cilist^2+1)^(N-1))) * (cilist^(N-3))
  cdfinverse <- 1-(cumsum(pdffunction)/sum(pdffunction)) # inverse of cdf
  criticalCI <- cilist[min(which(cdfinverse<=alpha))]  # find the threshold CI

  pval <- 0
  indices <- which(cilist>=CI)  # find the CI values larger than our CI
  if (length(indices)>0){pval <- cdfinverse[indices[1]]} # estimate the p-value

  output <- data.frame(CI,N,criticalCI,pval)
  return(output)}
