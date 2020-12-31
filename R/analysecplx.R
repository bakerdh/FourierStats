#' analysecplx: wrapper function that decides which statistical test to run for a given data set
#' the expected input is a data frame in long format:
#' one column should contain complex valued data, or two columns headed x and y contain the real and imaginary parts
#' if a column is headed 'group' (or any word beginning with a g) is included, this will be interpreted as level data for the independent variable
#' if a column is headed 'participant' (or any word beginning with a p) is included, this will be interpreted as a random factor
#'
#' the function first runs the condition index test for each level of the IV
#' if the condition index test is non-significant for all levels, a T-squared-circ test is run for one- and two-sample designs, or an ANOVA-squared-circ test otherwise
#' if the condition index test is significant for any level, a T-squared or MANOVA test is run instead
#' the paired-samples Mahalanobis effect size statistic is also calculated
#' @export
analysecplx <- function(input){

output <- NULL

# first check inputs and get into sensible format
quiterror <- 0
colheadings <- tolower(colnames(input))
# search for columns starting with x, y, g and p
xcol <- 0
ycol <- 0
groupcol <- 0
randcol <- 0
for (n in 1:ncol(input)){
  temp <- substr(colheadings[n],start=1,stop=1)
  if (pmatch('x',temp,nomatch=0)==1){xcol <- n}
  if (pmatch('y',temp,nomatch=0)==1){ycol <- n}
  if (pmatch('g',temp,nomatch=0)==1){groupcol <- n}
  if (pmatch('p',temp,nomatch=0)==1){randcol <- n}
}

# look for a column containing complex data
findcplx <- 0*ncol(input)
for (n in 1:ncol(input)){if (is.complex(input[1,n])){findcplx[n] <- 1}}
if (sum(findcplx)>0){cplxdata <- input[,which(findcplx==1)[1]]}

# if no complex data found, look for columns headed x and y
if (sum(findcplx)==0){
  if (xcol>0 && ycol>0){cplxdata <- complex(real=input[,xcol], imaginary=input[,ycol])}
  if (xcol==0 && ycol==0){quiterror <- 1
  errormsg <- 'No data found! Please include either a column of complex values, or two columns headed x and y'}
}

if (quiterror==0){
# determine study design, number of levels, within or between subjects
if (groupcol>0){grouping <- as.factor(input[,groupcol])}
if (groupcol==0){grouping <- as.factor(rep(1,nrow(input)))}
ngroups <- nlevels(grouping)
isRM <- 0
if (randcol>0){isRM <- 1
participants <- as.factor(input[,randcol])}

# do condition index test for each unique condition
grouplevs <- unique(grouping)
CIresults <- matrix(0,nrow=ngroups,ncol=4)
for (n in 1:ngroups){
  i <- which(grouping==grouplevs[n])
  CIresults[n,] <- unlist(CI.test(cplxdata[i]))
}
# compare condition index test for all groups to Bonferroni corrected alpha level
assumptionsmet <- 1
if (min(CIresults[,4])<(0.05/ngroups)){assumptionsmet <- 0}

# run either T-squared, T-squared-circ, MANOVA or ANOVA-squared-circ

data <- data.frame(Re(cplxdata),Im(cplxdata))

# run T-squared-circ
if (assumptionsmet==1 && ngroups<3){
  if (ngroups==1){results <- tsqc.test(data)}
  if (ngroups>1){
    dataA <- data[which(grouping==grouplevs[1]),]
    dataB <- data[which(grouping==grouplevs[2]),]
    results <- tsqc.test(dataA,dataB,paired=isRM)}
}

# run T-squared
if (assumptionsmet==0 && ngroups<3){
  if (ngroups==1){results <- tsq1.test(data)}
  if (ngroups>1){
    dataA <- data[which(grouping==grouplevs[1]),]
    dataB <- data[which(grouping==grouplevs[2]),]
  if (isRM==0){
    results <- Hotelling::hotelling.test(x=dataA,y=dataB)}
  if (isRM==1){
    results <- tsq1.test(dataA-dataB)   # take the difference between groups and run a 1 sample test
  }
}
}

# run ANOVA-squared-circ
if (assumptionsmet==1 && ngroups>2){

  if (isRM==0){
    print('Running independent ANOVA circ')
    dataAOV <- data.frame(cplxdata,grouping)
    colnames(dataAOV) <- c('simdata','grouplabels')
    results <- anovacirc.test(dataAOV)
  }
  if (isRM==1){
    print('Running paired samples ANOVA circ')
    dataAOV <- data.frame(cplxdata,grouping,participant)
    colnames(dataAOV) <- c('simdata','grouplabels','participant')
    results <- anovacircR.test(dataAOV)
  }
}

# run MANOVA
if (assumptionsmet==0 && ngroups>2){

  if (isRM==0){
    print('Running independent MANOVA')
  dataforManova <- data.frame(grouping,data[,1],data[,2])
  colnames(dataforManova) <- c("Group","Real","Imaginary")
  outcome <- cbind(dataforManova$Real, dataforManova$Imaginary)
  results <- summary(manova(outcome ~ Group, data = dataforManova))
  }
  if (isRM==1){
    print('Running repeated measures MANOVA')
    dataforManova <- data.frame(grouping,participant,data[,1],data[,2])
    colnames(dataforManova) <- c("Group","Participant","Real","Imaginary")
    outcome <- cbind(dataforManova$Real, dataforManova$Imaginary)
    # manmod <- lm(outcome ~ Group, data=dataforManova)
    results <- summary(manova(outcome ~ Group + Error(1/Participant), data = dataforManova))
    # manmod <- lm(outcome ~ Group + (1/Participant), data=dataforManova)
    # results <- summary(car::Manova(manmod, data = dataforManova))
  }

}

print(results)

# calculate pairwise Mahalanobis distance between each pair of conditions


# collate output



}
if (quiterror==1){output <- errormsg}

return(output)}
