Documentation for FourierStats Matlab package

The FourierStats package implements several statistical tests for use with bivariate data such as the complex (real and imaginary) components of the Fourier spectrum. Included functions implement the T-squared-circ test of Victor & Mast (1991), and the condition index and ANOVA-squared-circ tests described by Baker (2021). You should add the enclosing folder to your Matlab path to have access to all of the functions. Each function has documentation on the first few lines of the file, which can be accessed using the help function.

Functions include:

 - analysecplx: a function to analyse complex Fourier data, following the guidelines from the Baker (2021) paper

 - tsqc_test: implementation of the T-squared-circ test of Victor & Mast (1991)

 - tsqh_test: Hotelling's T-squared test (one and two sample versions)

 - CI_test: condition index test, to test the assumptions of the T-squared-circ test

 - anovacirc_test: implementation of one-way between subjects and repeated measures ANOVA-squared-circ test

 - amperrors: function to calculate error bars for amplitudes incorporating coherent averaging

 - pairwisemahal: calculates the pairwise Mahalanobis distance between two group means

 - clustercorrect: implements cluster correction described by Maris & Oostenveld (2007) for multivariate T-tests

 - getel: helper function that calculates the bounding ellipse for a cloud of points

 - fftshift: helper function that performs the quadrant shift of a 2D Fourier spectrum

Example data sets are contained in the files Hwangdata.mat and humanSSVEPdata.mat.

Package available from: https://github.com/bakerdh/FourierStats
For further details see: http://arxiv.org/abs/2101.04408
DHB 12/01/21
