# Adaptive trial designs with selection of composite endpoints and sample size reassessment

This repository contains the R package and main code to reproduce the results in ``Adaptive clinical trial designs with selection of composite endpoints and sample size reassessment’’ by M. Bofill Roig,  G. Gómez, M. Posch and F. Koenig.

In this work, we consider two-arm randomized controlled trials with a primary composite binary endpoint and an endpoint that consists only of the clinically most important component of the composite endpoint. We propose a trial design that allows an adaptive modification of the primary endpoint based on blinded information obtained at an interim analysis.


## R-package **eselect**

You can install the development version of the R package from GitHub:

``` r
# install.packages("devtools")
library(devtools)
install_github("MartaBofillRoig/eselect")
```

The R package contains the following functions:

- eselect: Endpoint selection and sample size reassessment for composite endpoints based on blinded data

- eselect_ub: Endpoint selection and sample size reassessment for composite endpoints based on unblinded data

- eselectsim: Simulation trials with endpoint selection and sample size reassessment for composite endpoints based on blinded data

- eselectsim_ub: Simulation trials with endpoint selection and sample size reassessment for composite endpoints based on unblinded data


## Reproducible research

This repository also contains the source files of the preprint (XXX). Specifically, 

- In the folder CODE_paper/example, there is the source code for reproduce the illustration; 
- In the folder CODE_paper/simulations, there is the code to reproduce the simulation study in the paper.
