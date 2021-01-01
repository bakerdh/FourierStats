#' anovacirc.test: two-dimensional analysis of variance using complex data
#' an extension of the logic of the T-squared circ statistic of Victor & Mast (1991)
#' this is a one-way implementation of the test - if participant IDs are supplied a repeated measures version is conducted
#'
#' Inputs--
#' data:  this can either be a vector of complex numbers, or a matrix
#'        if it is a matrix, either the first column contains complex values, or the first two columns are the x and y (real and imaginary) values of the DV
#'        if there are further columns, these are treated as the group labels and participant IDs
#' group: condition labels indicating the level of the independent variable that each data point corresponds to
#'        this is an optional input, but if it is supplied it supercedes values from the matrix
#' participant: variable storing participant/subject IDs (i.e. the random factor) for repeated measures analysis
#'        this is an optional input, but if it is supplied it supercedes values from the matrix
#'
#' see Baker (2021) for further details
#' @export
anovacirc.test <- function(data, group=NULL, participant=NULL){

  grouplabels <- NULL
  participantlabels <- NULL

  # if the input data is a vector, it must be complex Fourier components
  if (is.null(ncol(data))){datavals <- data}

  # if the input is not a vector, then work out what it contains based on its dimensions
  if (!is.null(ncol(data))){d <- dim(data)

  if (is.complex(data[,1])){datavals <- data[,1]
  if (d[2]>1){grouplabels <- data[,2]}
  if (d[2]>2){participantlabels <- data[,3]}}

  if (!is.complex(data[,1])){datavals <- complex(real=data[,1],imaginary=data[,2])
  if (d[2]>2){grouplabels <- data[,3]}
  if (d[2]>3){participantlabels <- data[,4]}}
  }

  if (!is.null(group)){grouplabels <- group}
  if (!is.null(participant)){participantlabels <- participant}

  grouplabels <- as.factor(grouplabels)
  factorlist <- levels(grouplabels)

  grandmean <- mean(datavals)


  # if no participant labels have been supplied, run a between-subjects ANOVA-circ
  if (is.null(participantlabels)){
  groupmeans <- NULL
  for (n in 1:nlevels(grouplabels)){groupmeans[n] <- mean(datavals[which(grouplabels==factorlist[n])])}

  SSM <- 0
  for (n in 1:nlevels(grouplabels)){SSM <- SSM + length(which(grouplabels==factorlist[n]))*abs(groupmeans[n]-grandmean)^2}

  dfM <- 2*(nlevels(grouplabels)-1)

  allresiduals <- 1:length(datavals)
  for (n in 1:nlevels(grouplabels)){allresiduals[which(grouplabels==factorlist[n])] <- abs(datavals[which(grouplabels==factorlist[n])] - groupmeans[n])}

  SSR <- sum(allresiduals^2)
  dfR <- 2*(length(datavals)-nlevels(grouplabels))
  MSM <- SSM/dfM
  MSR <- SSR/dfR
  Fratio <- MSM/MSR
  p.value <- pf(Fratio,df1=dfM,df2=dfR,lower.tail=FALSE)

  output <- data.frame(Fratio,p.value,SSM,SSR,dfM,dfR,MSM,MSR)
}

  # if participant labels have been supplied, run a repeated measures ANOVA-circ
  if (!is.null(participantlabels)){
  participantlabels <- as.factor(participantlabels)
  participantlist <- levels(participantlabels)

  SSW <- 0
  for (n in 1:nlevels(participantlabels)){SSW <- SSW + sum(abs(datavals[which(participantlabels==participantlist[n])] - mean(datavals[which(participantlabels==participantlist[n])]))^2)}
  dfW <- 2*(nlevels(participantlabels)*(nlevels(grouplabels)-1))

  groupmeans <- NULL
  for (n in 1:nlevels(grouplabels)){groupmeans[n] <- mean(datavals[which(grouplabels==factorlist[n])])}
  SSM <- 0
  for (n in 1:nlevels(grouplabels)){SSM <- SSM + length(which(grouplabels==factorlist[n]))*abs(groupmeans[n]-grandmean)^2}
  dfM <- 2*(nlevels(grouplabels)-1)

  SSR <- SSW - SSM
  dfR <- dfW - dfM

  MSM <- SSM/dfM
  MSR <- SSR/dfR
  Fratio <- MSM/MSR
  p.value <- pf(Fratio,df1=dfM,df2=dfR,lower.tail=FALSE)

  output <- data.frame(Fratio,p.value,SSW,SSM,SSR,dfW,dfM,dfR,MSM,MSR)
  }

  return(output)
}
