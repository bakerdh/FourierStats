#' clustercorrect: implements non-parametric cluster correction method described by Maris & Oostenveld (2007) J Neurosci Methods, 164: 177-190, doi: 10.1016/j.jneumeth.2007.03.024
#' this version can deal with complex data using the T-squared or T-squared-circ statistics
#'
#' The input datax should be an N (participants) x m (sensors/timepoints/locations) matrix
#' The optional input datay can either be a matrix of the same dimensions, or N2 x m (for unpaired designs with unbalanced samples), or a single value to compare to datax
#'
#' Legal values of the testtype input are 1: t-test, 2: Hotelling's T-squared, 3: Victor & Mast's T-squared-circ
#' Complex values are required for options 2 & 3. If complex values are passed in for test 1 (t-test), they will be converted to amplitudes
#'
#' if the measures are from adjacent time points, an adjacency matrix will be automatically constructed
#' if the measures are from different spatial locations (e.g. an electrode montage) you should calculate an adjacency matrix yourself
#' this should be an m x m matrix (where m is the number of sensors), where values are zero for non-adjacent sensor pairs, and 1 for adjacent pairs
#'
#' outputs an object containing the indices of all significant clusters
#' @export
clustercorrect <- function(datax,datay=NULL,adjacencymatrix=NULL,testtype=1,paired=TRUE,clustformthresh=0.05,clustthresh=0.05,nresamples=10000) {
  clustout <- NULL
  sigclustercount <- 0
  clusterpout <- NULL
  clusterpoints <- list()
  nulldist <- NULL
  N <- nrow(datax)  # number of observations (e.g. participants)
  m <- ncol(datax)  # number of measurements (e.g. sensors, voxels, timepoints)

  # if complex data have been passed in and t-test selected, convert to amplitudes
  if (is.complex(datax) && testtype==1){datax <- abs(datax)}
  if (is.complex(datay) && testtype==1){datay <- abs(datay)}

  # if two data sets have been passed in, note the size of the second matrix
  if (!pracma::isempty(datay)){N2 <- nrow(datay)}

  # if it's a paired design, subtract the two data sets and delete the second one
  if (paired==TRUE && !pracma::isempty(datay)){
    datax <- datax - datay
    datay <- NULL}

  if (length(adjacencymatrix)==0){
    # if no adjacency matrix has been passed in, we assume that subsequent observations are adjacent
    adjacencymatrix <- matrix(0,nrow=m,ncol=m)
    for (n in 1:(m-1)){adjacencymatrix[n,n+1] <- 1}
    for (n in 2:m){adjacencymatrix[n,n-1] <- 1}
  }

  allp <- NULL
  allt <- NULL
  # loop through all elements (sensors/time points) and calculate the test statistic and p-value
  for (n in 1:m){

    if (pracma::isempty(datay)){
      # one-sample or paired test
      if (testtype==1){
        output <- t.test(datax[,n])  # univariate t-test
        allp[n] <- output$p.value
        allt[n] <- output$statistic}
      if (testtype==2){
        output <- tsqh.test(datax[,n])  # one-sample Hotelling's t-squared
        allp[n] <- output$pval
        allt[n] <- output$tsq}
      if (testtype==3){
        output <- tsqc.test(datax[,n])  # T-squared-circ
        allp[n] <- output$pval
        allt[n] <- output$tsqc}
    }

    # independent samples test
    if (!pracma::isempty(datay) && is.matrix(datay)){
      if (testtype==1){
        output <- t.test(datax[,n],datay[,n]) # independent t-test
        allp[n] <- output$p.value
        allt[n] <- output$statistic}
      if (testtype==2){
        output <- tsqh.test(datax[,n],datay[,n],paired=FALSE)  # independent T-squared
        allp[n] <- output$pval
        allt[n] <- output$tsq}
      if (testtype==3){
        output <- tsqc.test(datax[,n],datay[,n],paired=FALSE)  # independent T-squared-circ
        allp[n] <- output$pval
        allt[n] <- output$tsqc}
    }
  }

  allp[which(is.na(allp))] <- 1

  # now generate a list of clusters of adjacent significant elements
  hvect <- rep(0,m)
  hvect[which(allp<clustformthresh)] <- 1
  hmat <- matrix(hvect,nrow=m,ncol=m)
  hmat2 <- hmat * t(hmat)
  clustprod <- hmat2 * adjacencymatrix

  clusterlist <- list()
  nclusts <- 0
  for (n in 1:m){
    temp = sum(clustprod[n,])
    if (temp>0){
      i <- which(clustprod[n,]>0)
      tempclust <- c(n, i)
      clustprod[n,] <- 0
      clustprod[,n] <- 0
      cnt <- 0
      while (cnt<length(tempclust)){
        cnt <- cnt + 1
        temp = sum(clustprod[tempclust[cnt],])
        if (temp>0){
          i <- which(clustprod[tempclust[cnt],]>0)
          tempclust <- c(tempclust, i)
          clustprod[tempclust[cnt],] <- 0
          clustprod[,tempclust[cnt]] <- 0
        }
      }
      nclusts <- nclusts + 1
      clusterlist[[nclusts]] <- sort(unique(tempclust))
    }
  }

  if (nclusts>0){
    sumtvals <- NULL
    for (cc in 1:nclusts){
      Cindices <- as.numeric(clusterlist[[cc]])
      sumtvals[cc] <- sum(allt[Cindices])}

    i <- which(sumtvals==max(sumtvals))  # find the largest cluster
    maxcluster <- clusterlist[[i]]  # store the largest cluster

    if (length(maxcluster)>1){
    # build a null distribution by permuting signs/group labels
    for (n in 1:nresamples){
      tsum <- 0
      if (pracma::isempty(datay)){randsigns <- (round(runif(N))*2)-1}
      if (!pracma::isempty(datay) && is.matrix(datay)){randgroups <- permute(c(1:N,(N+1):(N+N2)))}
      for (elcounter in 1:length(maxcluster)){

        if (pracma::isempty(datay)){  # one-sample or paired test
          tempdataA <- datax[,maxcluster[elcounter]]*randsigns
          if (testtype==1){
            output <- t.test(tempdataA)
            tsum <- tsum + output$statistic}
          if (testtype==2){
            output <- tsqh.test(tempdataA)
            tsum <- tsum + output$tsq}
          if (testtype==3){
            output <- tsqc.test(tempdataA)
            tsum <- tsum + output$tsqc}
        }

        if (!pracma::isempty(datay) && is.matrix(datay)){  # independent samples test
          tempdata <- c(datax[,maxcluster[elcounter]],datay[,maxcluster[elcounter]])
          tempdataA <- tempdata[randgroups(1:N)]
          tempdataB <- tempdata[randgroups((N+1):(N+N2))]
          if (testtype==1){
            output <- t.test(tempdataA,tempdataB)
            tsum <- tsum + output$statistic}
          if (testtype==2){
            output <- tsqh.test(tempdataA,tempdataB,paired=FALSE)
            tsum <- tsum + output$tsq}
          if (testtype==3){
            output <- tsqc.test(tempdataA,tempdataB,paired=FALSE)
            tsum <- tsum + output$tsqc}
        }
      }
      nulldist[n] <- tsum
    }
    }

    # compare each cluster to the null distribution, retain the significant ones
    clusterps <- 0*(1:nclusts)
    if (length(nulldist)==nresamples){
    for (cc in 1:nclusts){
      i <- which(abs(nulldist)>abs(sumtvals[cc]))
      if (!isempty(i)){clusterps[cc] <- length(i)/nresamples}
      if (clusterps[cc]<clustthresh){
        sigclustercount <- sigclustercount + 1
        clusterpoints[[sigclustercount]] <- clusterlist[[cc]]
        clusterpout[sigclustercount] <- clusterps[cc]
        }
    }}
  }
clustout <- NULL
clustout$clusterpoints <- clusterpoints
clustout$nclusters <- sigclustercount
clustout$pvals <- clusterpout

return(clustout)
}
