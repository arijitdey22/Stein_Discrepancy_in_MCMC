# Stein_Discrepancy_in_MCMC

This repository contains all the Rcodes, .Rdata file and the required plots in .pdf format that has been used in the prject. Below we describe the files as they have been created:

* 00_Draw_MVN_MCMC.R: The .R file which is sourced in almost every other parts of code and contains code for drawing MCMC samples that has a multivariate normal distribution as the stationary distribution.
2. 00_Stein_Thin_Indices.cpp: The C++ file that contains code for Stein-thinning of a MCMC samples with the target distribution being the standard normal distribution of appropriate order.
3. 01_KSD_MVN.pdf: The .pdf file that contains the plot of KSD and computation time with varrying values of sample size.
4. 01_KSD_MVN.R: The .R file that contains code for calculating KSDs and computation times mentioned in 01_KSD_MVN.pdf.
5. 01_KSD_MVN.Rdata: The .Rdata file that contains information to plot the plot from 01_KSD_MVN.pdf and obtained from 01_KSD_MVN.R.
6. 02_KSD_GM_MVN.pdf 
