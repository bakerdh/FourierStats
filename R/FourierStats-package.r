#' Documentation for FourierStats R package
#' The FourierStats package implements several statistical tests for use with bivariate data such as the complex (real and imaginary) components of the Fourier spectrum. Included functions implement the T-squared-circ test of Victor & Mast (1991), and the condition index and ANOVA-squared-circ tests described by Baker (2021).
#' Functions include:
#'  - analysecplx: a function to analyse complex Fourier data, following the guidelines from the Baker (2021) paper
#'  - tsqc.test: implementation of the T-squared-circ test of Victor & Mast (1991)
#'  - tsqh.test: Hotelling's T-squared test (one and two sample versions)
#'  - CI.test: condition index test, to test the assumptions of the T-squared-circ test
#'  - anovacirc.test: implementation of one-way between subjects and repeated measures ANOVA-squared-circ test
#'  - amperrors: function to calculate error bars for amplitudes incorporating coherent averaging
#'  - pairwisemahal: calculates the pairwise Mahalanobis distance between two group means
#'  - clustercorrect: implements cluster correction described by Maris & Oostenveld (2007) for multivariate T-tests
#'  - getel: helper function that calculates the bounding ellipse for a cloud of points
#'  - fftshift: helper function that performs the quadrant shift of a 2D Fourier spectrum
#'
#'  package available from: https://github.com/bakerdh/FourierStats
#'  for further details see: http://arxiv.org/abs/2101.04408
#'  DHB 12/01/21
