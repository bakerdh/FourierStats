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
#' the Mahalanobis distance effect size statistic (D) is also calculated
#' @export
analysecplx <- function(data, group=NULL, participant=NULL){

output <- NULL
complexdata <- NULL
grouplabels <- NULL
participantlabels <- NULL

# if the input data is a vector, it must be complex Fourier components
if (is.null(ncol(data))){complexdata <- data}

# if the input is not a vector, then work out what it contains based on its dimensions
if (!is.null(ncol(data))){d <- dim(data)

if (is.complex(data[,1])){complexdata <- data[,1]
if (d[2]>1){grouplabels <- data[,2]}
if (d[2]>2){participantlabels <- data[,3]}}

if (!is.complex(data[,1])){complexdata <- complex(real=data[,1],imaginary=data[,2])
if (d[2]>2){grouplabels <- data[,3]}
if (d[2]>3){participantlabels <- as.factor(data[,4])}}
}

if (!is.null(group)){grouplabels <- group}
if (!is.null(participant)){participantlabels <- as.factor(participant)}

if (is.null(grouplabels)){grouplabels <- rep(1,length(complexdata))}
isRM <- 1
if (is.null(participantlabels)){isRM <- 0}

grouplabels <- as.factor(grouplabels)
factorlist <- levels(grouplabels)
ngroups <- nlevels(grouplabels)
nobservations <- length(complexdata)

message(paste('Design has',ngroups,'levels, with',nobservations/ngroups,'observations per level'))
if (!is.null(participantlabels)){message('Design is repeated measures')}

# run the condition index test for each group
CIresults <- matrix(0,nrow=ngroups,ncol=4)
for (n in 1:ngroups){
  i <- which(grouplabels==factorlist[n])
  CIresults[n,] <- unlist(CI.test(complexdata[i]))
}
nsigCI <- length(which(CIresults[,4]<(0.05/ngroups)))  # compare condition index test for all groups to Bonferroni corrected alpha level
message(paste('The condition index was significant for',nsigCI,'out of',ngroups,'levels'))
assumptionsmet <- 1
if (nsigCI>0){assumptionsmet <- 0}
output$ngroups <- ngroups
output$nsigCI <- nsigCI



# run T-squared-circ
if (assumptionsmet==1 && ngroups<3){
  if (ngroups==1){
    message('Running a one-sample T-squared-circ test')
    output$testtype <- 'One-sample T-squared-circ test'
    results <- tsqc.test(complexdata)}
  if (ngroups>1){
    if (isRM==0){message('Running an independent T-squared-circ test')
      output$testtype <- 'Independent T-squared-circ test'}
    if (isRM==1){message('Running a repeated measures T-squared-circ test')
      output$testtype <- 'Repeated measures T-squared-circ test'}
    dataA <- complexdata[which(grouplabels==factorlist[1])]
    dataB <- complexdata[which(grouplabels==factorlist[2])]
    results <- tsqc.test(dataA,dataB,paired=isRM)}
  message(paste('The test statistic T^2_circ =',round(results$tsqc,digits=2)))
  message(paste('The equivalent F-ratio with',results$df1,'and',results$df2,'degrees of freedom is F =',round(results$Fratio,digits=2)))
  output$teststat <- results$tsqc
  output$Fratio <- results$Fratio
  output$df1 <- results$df1
  output$df2 <- results$df2
  output$p.value <- results$p.value
}

# run T-squared
if (assumptionsmet==0 && ngroups<3){
  if (ngroups==1){
    message('Running a one-sample T-squared test')
    output$testtype <- 'One-sample T-squared test'
    results <- tsq1.test(complexdata)
    message(paste('The test statistic T^2 =',round(results$tsq,digits=2)))
    message(paste('The equivalent F-ratio with',results$df1,'and',results$df2,'degrees of freedom is F =',round(results$Fratio,digits=2)))
    output$teststat <- results$tsq
    output$Fratio <- results$Fratio
    output$df1 <- results$df1
    output$df2 <- results$df2
    output$p.value <- results$p.value
    }
  if (ngroups>1){
    dataA <- complexdata[which(grouplabels==factorlist[1])]
    dataB <- complexdata[which(grouplabels==factorlist[2])]
    if (isRM==0){
      message('Running an independent T-squared test')
      output$testtype <- 'Independent T-squared test'
      results <- Hotelling::hotelling.test(x=data.frame(Re(dataA),Im(dataA)),y=data.frame(Re(dataB),Im(dataB)))
      message(paste('The test statistic T^2 =',round(results$stats$statistic,digits=2)))
      output$teststat <- results$stats$statistic
      output$Fratio <- results$stats$statistic * results$stats$m
      output$df1 <- results$stats$df[1]
      output$df2 <- results$stats$df[2]
      output$p.value <- results$pval
      message(paste('The equivalent F-ratio with',output$df1,'and',output$df2,'degrees of freedom is F =',round(output$Fratio,digits=2)))
    }
    if (isRM==1){
      message('Running a repeated measures T-squared test')
      output$testtype <- 'Repeated measures T-squared test'
      results <- tsq1.test(dataA-dataB)   # take the difference between groups and run a 1 sample test
      message(paste('The test statistic T^2 =',round(results$tsq,digits=2)))
      message(paste('The equivalent F-ratio with',results$df1,'and',results$df2,'degrees of freedom is F =',round(results$Fratio,digits=2)))
      output$teststat <- results$tsq
      output$Fratio <- results$Fratio
      output$df1 <- results$df1
      output$df2 <- results$df2
      output$p.value <- results$p.value
      }
  }
}

# run ANOVA-squared-circ
if (assumptionsmet==1 && ngroups>2){
  if (isRM==0){message('Running an independent ANOVA-squared-circ test')
    output$testtype <- 'Independent ANOVA-squared-circ test'}
  if (isRM==1){message('Running a repeated measures ANOVA-squared-circ test')
    output$testtype <- 'Repeated measures ANOVA-squared-circ test'}
    results <- anovacirc.test(complexdata,grouplabels,participantlabels)
    message(paste('The F-ratio with',results$dfM,'and',results$dfR,'degrees of freedom is F =',round(results$Fratio,digits=2)))
    output$Fratio <- results$Fratio
    output$df1 <- results$dfM
    output$df2 <- results$dfR
    output$p.value <- results$p.value
}

# run MANOVA
if (assumptionsmet==0 && ngroups>2){

  if (isRM==0){
    message('Running an independent MANOVA')
    output$testtype <- 'Independent MANOVA'
    dataforManova <- data.frame(grouplabels,Re(complexdata),Im(complexdata))
    colnames(dataforManova) <- c("Group","Real","Imaginary")
    results <- MANOVA.RM::MANOVA.wide(cbind(dataforManova$Real, dataforManova$Imaginary) ~ Group, data = dataforManova)
  }
  if (isRM==1){
    message('Running a repeated measures MANOVA')
    dataforManova <- data.frame(grouplabels,participantlabels,Re(complexdata),Im(complexdata))
    colnames(dataforManova) <- c("Group","Participant","Real","Imaginary")
    results <- MANOVA.RM::multRM(cbind(dataforManova$Real, dataforManova$Imaginary) ~ Group, data=dataforManova, subject="Participant", within="Group")
    output$testtype <- 'Repeated measures MANOVA'
    }
  output$teststat <- results$MATS
  output$p.value <- results$resampling[2]
}

if (output$p.value<0.05){message('The test was significant at p < 0.05')}
if (output$p.value>=0.05){message('The test was not significant at p < 0.05')}

if (ngroups==1){
  # calculate pointwise Mahalanobis distance relative to the origin
  temp <- data.frame(Re(complexdata),Im(complexdata))
  D <- sqrt(stats::mahalanobis(c(0,0), colMeans(temp), cov(temp)))
}
if (ngroups > 1){
  # calculate pairwise Mahalanobis distance between each pair of conditions
  mahala_sq <- HDMD::pairwise.mahalanobis(data.frame(Re(complexdata),Im(complexdata)),grouplabels)
  D <- max(sqrt(mahala_sq$distance))  # choose the largest effect size
}
message(paste('The effect size (Mahalanobis distance) is D =',round(D,digits=2)))

output$D <- D

return(output)}



