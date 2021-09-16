# FourierStats
Statistical Tests for Periodic Fourier Data from Neuroscience Experiments

This github repository contains a package called FourierStats, implemented in R and Matlab. The R version can be installed with the following command:

devtools::install_github("bakerdh/FourierStats")

and activated with:

library(FourierStats)

This requires the devtools package to be installed.

All scripts are stored in the 'R' subdirectory, and there is a vignette demonstrating function usage in the doc directory. Also included in the 'manuscript' subdirectory is the Rmarkdown file (manuscript.Rmd) used to generate the paper "Statistical analysis of periodic data in neuroscience", as well as .tex and .pdf versions, and all figures and other resources.

The Matlab toolbox is stored in the FourierStatsMatlab directory. This should be downloaded and added to the Matlab path either using the Set Path GUI, or the addpath function.

The packages are available from: https://github.com/bakerdh/FourierStats

The full citation for the study is:

Baker, D.H. (2021). Statistical analysis of periodic data in neuroscience. Neurons, Behavior, Data Analysis and Theory, 5(3): 1-18. http://dx.doi.org/10.51628/001c.27680

