#' anovacirc.test: two-dimensional analysis of variance using complex data
#' an extension of the logic of the t-squared circ statistic
#' this is the one-way, between subjects implementation of the test
#' the expected input is a long-format data frame with the following columns:
#' simdata: the data stored as complex values
#' grouplabels: condition labels indicating the level of the independent variable each data point corresponds to
#' see Baker (2021) for further details
#' @export
anovacirc.test <- function(data){

  grouplabels <- data$grouplabels
  datavals <- data$simdata
  factorlist <- levels(grouplabels)

  grandmean <- mean(datavals)

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
  pvalue <- pf(Fratio,df1=dfM,df2=dfR,lower.tail=FALSE)

  output <- NULL
  output$Fratio <- Fratio
  output$pvalue <- pvalue
  output$SSM <- SSM
  output$SSR <- SSR
  output$dfM <- dfM
  output$dfR <- dfR
  output$MSM <- MSM
  output$MSR <- MSR

  return(output)
}
