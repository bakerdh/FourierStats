#' anovacircR.test: two-dimensional analysis of variance using complex data
#' an extension of the logic of the t-squared circ statistic
#' this is the one-way, repeated measures implementation of the test
#' the expected input is a long-format data frame with the following columns:
#' simdata: the data stored as complex values
#' grouplabels: condition labels indicating the level of the independent variable each data point corresponds to
#' participant: participant ID number
#' see Baker (2021) for further details
#' @export
anovacircR.test <- function(data){

  grouplabels <- data$grouplabels
  factorlist <- levels(grouplabels)
  participantlabels <- data$participant
  participantlist <- levels(participantlabels)
  datavals <- data$simdata

  grandmean <- mean(datavals)

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
  pvalue <- pf(Fratio,df1=dfM,df2=dfR,lower.tail=FALSE)

  output <- NULL
  output$Fratio <- Fratio
  output$pvalue <- pvalue
  output$SSW <- SSW
  output$SSM <- SSM
  output$SSR <- SSR
  output$dfW <- dfW
  output$dfM <- dfM
  output$dfR <- dfR
  output$MSM <- MSM
  output$MSR <- MSR

  return(output)
}
