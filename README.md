# miaSim - Microbiome data simulation package

This package is aimed to provide simulations for analyses of microbiome data by using community models.
The main class for working with microbiome data in this package is `SummarizedExperiment`.

The latest version of the package includes the following simulations :

- The Hubbell Neutral Model Community Simulation: `miaSim::simulateHubbell()`
- The generalized Lotka-Volterra Model Community Time Series Simulation: `miaSim::simulateGLV()`
- The Self-Organised Instability Model Community Simulation: `miaSim::simulateSOI()`

## Contribution

Feel free to contribute.

## Installation
 
The package can be directly installed from R command line.

```{R}
devtools::install_github("microbiome/miaSim")
library(miaSim)
```
## References
G. (2020–2021). GitHub - gheysenemma/microsimR: A Simulator Package with a collection of basic community models 
that outputs count data of microbial communities. GitHub. https://github.com/gheysenemma/microsimR

## Acknowledgements

- Author: Yagmur Simsek

Under supervision of:
- Promotor and Mentor: Prof. Leo Lahti(University of Turku)
