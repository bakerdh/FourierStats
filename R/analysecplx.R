#' analysecplx: wrapper function that automatically decides which statistical test to run for a given data set
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
#' the function first runs the condition index test for each level of the IV
#' if the condition index test is non-significant for all levels, a T-squared-circ test is run for one- and two-sample designs, or an ANOVA-squared-circ test otherwise
#' if the condition index test is significant for any level, a T-squared test or MANOVA is run instead
#' if participant labels are included, a repeated measures test is conducted
#' the paired-samples Mahalanobis effect size statistic is also calculated
#' @export
analysecplx <- function(data, group=NULL, participant=NULL){

output <- NULL
complexdata <- NULL
grouplabels <- NULL
participantlabels <- NULL

if (!is.matrix(data)){complexdata <- data}

if (is.matrix(data)){d <- dim(data)

if (is.complex(data[,1])){complexdata <- data[,1]
if (d[2]>1){grouplabels <- data[,2]}
if (d[2]>2){participantlabels <- data[,3]}}

if (!is.complex(data[,1])){complexdata <- complex(real=data[,1],imaginary=data[,2])
if (d[2]>2){grouplabels <- data[,3]}
if (d[2]>3){participantlabels <- data[,4]}}
}

if (!is.null(group)){grouplabels <- group}
if (!is.null(participant)){participantlabels <- participant}

if (is.null(grouplabels)){grouplabels <- rep(1,length(complexdata))}
isRM <- 1
if (is.null(participantlabels)){isRM <- 0}

grouplabels <- as.factor(grouplabels)
factorlist <- levels(grouplabels)
ngroups <- nlevels(grouplabels)

# run the condition index test for each group
CIresults <- matrix(0,nrow=ngroups,ncol=4)
for (n in 1:ngroups){
  i <- which(grouplabels==factorlist[n])
  CIresults[n,] <- unlist(CI.test(complexdata[i]))
}
# compare condition index test for all groups to Bonferroni corrected alpha level
assumptionsmet <- 1
if (min(CIresults[,4])<(0.05/ngroups)){assumptionsmet <- 0}


# run T-squared-circ
if (assumptionsmet==1 && ngroups<3){
  if (ngroups==1){results <- tsqc.test(complexdata)}
  if (ngroups>1){
    dataA <- complexdata[which(grouplabels==factorlist[1])]
    dataB <- complexdata[which(grouplabels==factorlist[2])]
    results <- tsqc.test(dataA,dataB,paired=isRM)}
}

# run T-squared
if (assumptionsmet==0 && ngroups<3){
  if (ngroups==1){results <- tsq1.test(complexdata)}
  if (ngroups>1){
    dataA <- complexdata[which(grouplabels==factorlist[1])]
    dataB <- complexdata[which(grouplabels==factorlist[2])]
    if (isRM==0){
      results <- Hotelling::hotelling.test(x=data.frame(Re(dataA),Im(dataA)),y=data.frame(Re(dataB),Im(dataB)))}
    if (isRM==1){
      results <- tsq1.test(dataA-dataB)   # take the difference between groups and run a 1 sample test
    }
  }
}

# run ANOVA-squared-circ
if (assumptionsmet==1 && ngroups>2){
    results <- anovacirc.test(complexdata,grouplabels,participantlabels)
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









# print(results)


if (ngroups==1){
  # calculate pointwise Mahalanobis distance relative to the origin
  temp <- data.frame(Re(complexdata),Im(complexdata))
  D <- sqrt(stats::mahalanobis(c(0,0), colMeans(temp), cov(temp)))
}
if (ngroups > 1){
  # calculate pairwise Mahalanobis distance between each pair of conditions
  mahala_sq <- HDMD::pairwise.mahalanobis(data.frame(Re(complexdata),Im(complexdata)),grouplabels)
  D <- sqrt(mahala_sq$distance)
}

# collate output




# if (quiterror==1){output <- errormsg}

return(output)}



#
#
# # first check inputs and get into sensible format
# quiterror <- 0
# colheadings <- tolower(colnames(input))
# # search for columns starting with x, y, g and p
# xcol <- 0
# ycol <- 0
# groupcol <- 0
# randcol <- 0
# for (n in 1:ncol(input)){
#   temp <- substr(colheadings[n],start=1,stop=1)
#   if (pmatch('x',temp,nomatch=0)==1){xcol <- n}
#   if (pmatch('y',temp,nomatch=0)==1){ycol <- n}
#   if (pmatch('g',temp,nomatch=0)==1){groupcol <- n}
#   if (pmatch('p',temp,nomatch=0)==1){randcol <- n}
# }
#
# # look for a column containing complex data
# findcplx <- 0*ncol(input)
# for (n in 1:ncol(input)){if (is.complex(input[1,n])){findcplx[n] <- 1}}
# if (sum(findcplx)>0){cplxdata <- input[,which(findcplx==1)[1]]}
#
# # if no complex data found, look for columns headed x and y
# if (sum(findcplx)==0){
#   if (xcol>0 && ycol>0){cplxdata <- complex(real=input[,xcol], imaginary=input[,ycol])}
#   if (xcol==0 && ycol==0){quiterror <- 1
#   errormsg <- 'No data found! Please include either a column of complex values, or two columns headed x and y'}
# }
#
# if (quiterror==0){
# # determine study design, number of levels, within or between subjects
# if (groupcol>0){grouping <- as.factor(input[,groupcol])}
# if (groupcol==0){grouping <- as.factor(rep(1,nrow(input)))}
# ngroups <- nlevels(grouping)
# isRM <- 0
# if (randcol>0){isRM <- 1
# participants <- as.factor(input[,randcol])}
#
# # do condition index test for each unique condition
# grouplevs <- unique(grouping)
# CIresults <- matrix(0,nrow=ngroups,ncol=4)
# for (n in 1:ngroups){
#   i <- which(grouping==grouplevs[n])
#   CIresults[n,] <- unlist(CI.test(cplxdata[i]))
# }
# # compare condition index test for all groups to Bonferroni corrected alpha level
# assumptionsmet <- 1
# if (min(CIresults[,4])<(0.05/ngroups)){assumptionsmet <- 0}
#
# # run either T-squared, T-squared-circ, MANOVA or ANOVA-squared-circ
#
# data <- data.frame(Re(cplxdata),Im(cplxdata))
#
# # run T-squared-circ
# if (assumptionsmet==1 && ngroups<3){
#   if (ngroups==1){results <- tsqc.test(data)}
#   if (ngroups>1){
#     dataA <- data[which(grouping==grouplevs[1]),]
#     dataB <- data[which(grouping==grouplevs[2]),]
#     results <- tsqc.test(dataA,dataB,paired=isRM)}
# }
#
# # run T-squared
# if (assumptionsmet==0 && ngroups<3){
#   if (ngroups==1){results <- tsq1.test(data)}
#   if (ngroups>1){
#     dataA <- data[which(grouping==grouplevs[1]),]
#     dataB <- data[which(grouping==grouplevs[2]),]
#   if (isRM==0){
#     results <- Hotelling::hotelling.test(x=dataA,y=dataB)}
#   if (isRM==1){
#     results <- tsq1.test(dataA-dataB)   # take the difference between groups and run a 1 sample test
#   }
# }
# }
#
# # run ANOVA-squared-circ
# if (assumptionsmet==1 && ngroups>2){
#
#   if (isRM==0){
#     print('Running independent ANOVA circ')
#     dataAOV <- data.frame(cplxdata,grouping)
#     colnames(dataAOV) <- c('simdata','grouplabels')
#     results <- anovacirc.test(dataAOV)
#   }
#   if (isRM==1){
#     print('Running paired samples ANOVA circ')
#     dataAOV <- data.frame(cplxdata,grouping,participant)
#     colnames(dataAOV) <- c('simdata','grouplabels','participant')
#     results <- anovacircR.test(dataAOV)
#   }
# }
#
# # run MANOVA
# if (assumptionsmet==0 && ngroups>2){
#
#   if (isRM==0){
#     print('Running independent MANOVA')
#   dataforManova <- data.frame(grouping,data[,1],data[,2])
#   colnames(dataforManova) <- c("Group","Real","Imaginary")
#   outcome <- cbind(dataforManova$Real, dataforManova$Imaginary)
#   results <- summary(manova(outcome ~ Group, data = dataforManova))
#   }
#   if (isRM==1){
#     print('Running repeated measures MANOVA')
#     dataforManova <- data.frame(grouping,participant,data[,1],data[,2])
#     colnames(dataforManova) <- c("Group","Participant","Real","Imaginary")
#     outcome <- cbind(dataforManova$Real, dataforManova$Imaginary)
#     # manmod <- lm(outcome ~ Group, data=dataforManova)
#     results <- summary(manova(outcome ~ Group + Error(1/Participant), data = dataforManova))
#     # manmod <- lm(outcome ~ Group + (1/Participant), data=dataforManova)
#     # results <- summary(car::Manova(manmod, data = dataforManova))
#   }
#
# }
