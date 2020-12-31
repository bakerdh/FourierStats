#' Documentation for FourierStats R package
#' The FourierStats package implements several statistical tests for use with bivariate data such as the complex (real and imaginary) components of the Fourier spectrum. Included functions implement the T-squared-circ test of Victor & Mast (1991), and the condition index and ANOVA-squared-circ described by Baker (2021).
#' Functions include:
#'  - analysecplx: a function to analyse complex Fourier data, following the guidelines from the Baker (2021) paper
#'  - tsqc.test: implementation of the T-squared-circ test of Victor & Mast (1991)
#'  - tsq1.test: one-sample Hotelling's T-squared test
#'  - CI.test: condition index test, to test the assumptions of the T-squared-circ test
#'  - anovacirc.test: implementation of one-way between subjects ANOVA-squared-circ test
#'  - anovacircR.test: implementation of one-way repeated measures ANOVA-squared-circ test
#'  - amperrors: function to calculate error bars for amplitudes incorporating coherent averaging
#'  - getel: helper function that calculates the bounding ellipse for a cloud of points
#'  - fftshift: helper function that performs the quadrant shift of a 2D Fourier spectrum
#'
#'  DHB 26/12/20
