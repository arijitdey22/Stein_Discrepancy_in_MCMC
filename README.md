# Stein_discrepancy_in_MCMC

This repository contains all the Rcodes, .Rdata file and the required plots in .pdf format that have been used in the project. Below we describe the files as they have been created:

* **00_Draw_MVN_MCMC.R**: The .R file which is sourced in almost every other parts of code and contains code for drawing MCMC samples that has a multivariate normal distribution as the stationary distribution.
* **00_Stein_Thin_Indices.cpp**: The C++ file that contains code for Stein-thinning of a MCMC samples with the target distribution being the standard normal distribution of appropriate order.
* **01_KSD_MVN.pdf**: The .pdf file that is the plot of KSD and computation time with varrying values of sample size.
* **01_KSD_MVN.R**: The .R file that contains code for calculating KSDs and computation times mentioned in 01_KSD_MVN.pdf.
* **01_KSD_MVN.Rdata**: The .Rdata file that contains information to create the plot from 01_KSD_MVN.pdf and obtained from 01_KSD_MVN.R.
* **02_KSD_GM_MVN.pdf**: The .pdf file that is the plot of Stein-thinning of Multivariate normal samples obtained from MCMC.
* **02_KSD_GM_MVN.R**: The .R file that contains code for choosing the Stein-thinned indices.
* **03_KSD_vs_Stein_thinned_KSD.1.pdf**: The .pdf file that is the plot that compare the KSD of full sample with the KSD of the Stein-thinned sample.
* **03_KSD_vs_Stein_thinned_KSD.2.pdf**: The .pdf file that is the plot that compare the computation time to calcaulate the KSD of full and thinned sample and also the time required for Stein-thinning. 
* **03_KSD_vs_Stein_thinned_KSD.R**: The .R file that contains code for calculating the KSDs and the comoutation times needed for the above plots.
* **03_KSD_vs_Stein_thinned_KSD.Rdata**: The .Rdata file that contains information to create the plot from 01_KSD_MVN.pdf and obtained from 01_KSD_MVN.R.
* **04_Thin_comp.pdf**: The .pdf file that is the plot that compares the MSE of Full-sample and three types of thinned-sample estimate of the parameter of interest.
* **04_Thin_comp.R**: The .R file that contains code for computing the MSEs mentioned above.
* **04_Thin_comp.Rdata**: The .Rdata file that contains information to plot the create from 01_KSD_MVN.pdf and obtained from 01_KSD_MVN.R.
* **05_Thin_comp_Biased_MCMC.pdf**: The .pdf file that is the plot that compares the MSE of Full-sample and three types of thinned-sample estimate of the parameter of interest for a biased MCMC.
* **05_Thin_comp_Biased_MCMC.R**: The .R file that contains code for computing the MSEs mentioned above.
* **05_Thin_comp_Biased_MCMC.Rdata**: The .Rdata file that contains information to create the plot from 01_KSD_MVN.pdf and obtained from 01_KSD_MVN.R.
* **06_Thin_comp_Biased_MCMC_Examples.250.pdf**: The .pdf file that is the plot of densitites of four scenarios for biased MCMC with n = 250.
* **06_Thin_comp_Biased_MCMC_Examples.250.Rdata**: The .Rdata file that contains information to create the plot from 06_Thin_comp_Biased_MCMC_Examples.250.pdf and obtained from 06_Thin_comp_Biased_MCMC_Examples.R.
* **06_Thin_comp_Biased_MCMC_Examples.500.pdf**: The .pdf file that is the plot of densitites of four scenarios for biased MCMC with n = 500.
* **06_Thin_comp_Biased_MCMC_Examples.500.Rdata**: The .Rdata file that contains information to create the plot from 06_Thin_comp_Biased_MCMC_Examples.500.pdf and obtained from 06_Thin_comp_Biased_MCMC_Examples.R.
* **06_Thin_comp_Biased_MCMC_Examples.1000.pdf**: The .pdf file that is the plot of densitites of four scenarios for biased MCMC with n = 1000.
* **06_Thin_comp_Biased_MCMC_Examples.1000.Rdata**: The .Rdata file that contains information to create the plot from 06_Thin_comp_Biased_MCMC_Examples.1000.pdf and obtained from 06_Thin_comp_Biased_MCMC_Examples.R.
* **06_Thin_comp_Biased_MCMC_Examples.5000.pdf**: The .pdf file that is the plot of densitites of four scenarios for biased MCMC with n = 5000.
* **06_Thin_comp_Biased_MCMC_Examples.5000.Rdata**: The .Rdata file that contains information to create the plot from 06_Thin_comp_Biased_MCMC_Examples.5000.pdf and obtained from 06_Thin_comp_Biased_MCMC_Examples.R.
* **06_Thin_comp_Biased_MCMC_Examples.R**: The .R file that contains code for calculating the densities of the four scenarios mentioned in the above four plots.

