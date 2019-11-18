This repository contains software made available by [JanBiotech](https://janbiotech.com/) to estimate cell counts from limited dilution assays.

Simulations
-----------

The _simulation_ folder contains the R code for simulating such assays and plotting the results. Running these simulations trains intuition and provides the ground truth for testing of model implementations. The simulations currently implemented follow the protocol of the [QVOA assay](https://www.biorxiv.org/content/10.1101/018911v1) to estimate the latent reservoir of replication competent HIV. The scripts depend on the [data.table](https://github.com/Rdatatable/data.table/wiki) R package. Plotting depends on `ggplot2` and `showtext` packages. I use the Myriad Pro fonts, users can substitute their own in the `font_add()` function call in the `simPlots.Rnw` script.

Quantal dilution assay analyses
-------------------------

The R package for Bayesian analyses of quantal dilution assays (such as QVOA) is in the BayesQLD directory of this repository. It depends only on [Rcpp](http://www.rcpp.org/). It can be installed directly from GitHub by running `install_github("tonymugen/BayesLiDiCal/BayesQLD")`, which requires the `devtools` package.

`BayesQLD` is still in development, check back for newer versions.

