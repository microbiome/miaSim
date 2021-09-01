## miaSim - Microbiome data simulation package

<!-- badges: start -->

[![R-CMD-check-bioc-devel](https://github.com/microbiome/miaSim/workflows/R-CMD-check-bioc-devel/badge.svg)](https://github.com/microbiome/miaSim/actions)
[![R-CMD-check-bioc](https://github.com/microbiome/miaSim/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/microbiome/miaSim/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/microbiome/miaSim/branch/main/graph/badge.svg)](https://codecov.io/gh/microbiome/miaSim?branch=main)

<!-- badges: end -->

This package is aimed to provide simulations for analyses of microbiome data by using community models.
The main class for working with microbiome data in this package is `SummarizedExperiment`.

The latest version of the package includes the following simulations :

- The Hubbell Neutral Model community simulation: `miaSim::simulateHubbell()`
- The generalized Lotka-Volterra model community time series simulation: `miaSim::simulateGLV()`

# Contribution

Feel free to contribute.

## Installation
 
The package can be directly installed from R command line.

```{R}
devtools::install_github("microbiome/miaSim")
library(miaSim)
```

# Acknowledgements

- Author: Yagmur Simsek

Under supervision of:
- Promotor and Mentor: Prof. Leo Lahti(University of Turku)
