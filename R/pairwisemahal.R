#' pairwisemahal: calculates the pairwise (group to group) Mahalanobis distance between group means
#'
#' inputs:
#' data - can be either a vector of N complex numbers, a matrix of Nxm complex numbers, or an Nx2 matrix of real numbers
#'        if a vector of complex numbers or matrix of real numbers, the 'grouping' variable is also required
#'        if a matrix of complex numbers is supplied, the columns are interpreted as the groups
#' grouping - a variable containing group membership information for the data
#'
#' output:
#' an mxm data frame of Mahalanobis distances between the group means
#'     note that the Mahalanobis distance (D) is returned, not D^2
#' the row and column names of the data frame are the condition names from the grouping input variable (or are otherwise numeric)
#'
#' this function is part of the FourierStats package: https://github.com/bakerdh/FourierStats

#' @export
pairwisemahal <- function(data, grouping=NULL){

  d <- dim(data)
  if (!is.null(ncol(data))){if (d[1]<d[2]){data <- t(data)}}
  d <- dim(data)
  if (is.complex(data)){
    if (min(d)==1){       # it is a vector
      realdata <- data.frame(Re(data),Im(data))
    }
    if (min(d)>1){        # it is a matrix
      ngroups <- d(2)
      realdata <- matrix(0,nrow=d(1),ncol=2)
      grouping <- NULL
      for (n in 1:ngroups){
        realdata[((n-1)*d(1)+1):(n*d(1)),1] <- Re(data)
        realdata[((n-1)*d(1)+1):(n*d(1)),2] <- Im(data)
        grouping[((n-1)*d(1)+1):(n*d(1))] <- n
      }
    }
  }
  if (!is.complex(data)){realdata = data}

  grouplabels <- as.factor(grouping)
  factorlist <- levels(grouplabels)
  m = length(factorlist)

  centroids <- matrix(0,nrow=m,ncol=2)
  for (n in 1:m){centroids[n,] <- colMeans(realdata[which(grouplabels==factorlist[n]),])}

  distances <- matrix(0,nrow=m,ncol=m)
  for (n1 in 1:m){
    for (n2 in 1:m){
      CA <- cov(realdata[which(grouplabels==factorlist[n1]),])
      CB <- cov(realdata[which(grouplabels==factorlist[n2]),])
      N1 <- length(which(grouplabels==factorlist[n1]))
      N2 <- length(which(grouplabels==factorlist[n2]))
      pwcov <- ((N1-1)*CA + (N2-1)*CB)/(N1+N2-2)
      meandiff <- centroids[n1,] - centroids[n2,]
      distances[n1,n2] <- sqrt(t(meandiff) %*% solve(pwcov) %*% meandiff)
    }
  }

  output <- data.frame(distances)
  colnames(output) <- factorlist
  rownames(output) <- factorlist
  return(output)}
