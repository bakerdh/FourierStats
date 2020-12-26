#' helper function that calculates an ellipse bounding some data points
#' the outline of the ellipse is returned
#' @export
getel <- function(compdata){
  A <- cov(compdata)
  ctr    <- colMeans(compdata)
  RR     <- chol(A)                               # Cholesky decomposition
  angles <- seq(0, 2*pi, length.out=200)          # angles for ellipse
  ell    <- 1 * cbind(cos(angles), sin(angles)) %*% RR  # ellipse scaled with factor 1
  ellCtr <- sweep(ell, 2, ctr, "+")               # center ellipse to the data centroid
  eigVal  <- eigen(A)$values
  eigVec  <- eigen(A)$vectors
  eigScl  <- eigVec  %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
  xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
  yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
  ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles))
  ellRot  <- eigVec %*% t(ellBase)
  return(ellRot)
}
