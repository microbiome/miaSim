# Microbiome data simulation with miaSim

<!-- badges: start -->
[![R-CMD-check](https://github.com/microbiome/miaSim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/microbiome/miaSim/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check-bioc](https://github.com/microbiome/miaSim/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/microbiome/miaSim/actions/workflows/check-bioc.yml)
[![R-CMD-check-bioc-devel](https://github.com/microbiome/miaSim/actions/workflows/check-bioc-devel.yml/badge.svg)](https://github.com/microbiome/miaSim/actions/workflows/check-bioc-devel.yml)
[![Codecov test coverage](https://codecov.io/gh/microbiome/miaSim/branch/master/graph/badge.svg)](https://codecov.io/gh/microbiome/miaSim?branch=master)
[![Gitter](https://badges.gitter.im/microbiome/mia.svg)](https://gitter.im/microbiome/miaverse)

[![Watch on GitHub][github-watch-badge]][github-watch]
[![Star on GitHub][github-star-badge]][github-star]



<!-- badges: end -->

This R package can be used to simulate (longitudinal) data for the
benchmarking and analysis of quantitative models of microbial
communities. The package is part of
[miaverse](microbiome.github.io), and is based on the
`(Tree)SummarizedExperiment` data container.

The package [website](https://microbiome.github.io/miaSim/) provides
further tutorials and references for the [implemented
models](https://microbiome.github.io/miaSim/reference/index.html),
which include Hubbell's neutral model, the generalized Lotka-Volterra
model, and the self-organised instability model.


## Installation

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("miaSim")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("miaSim")
```

### miaSimShiny

miaSim is accompanied by [miaSimShiny](https://github.com/gaoyu19920914/miaSimShiny), which allows users to explore the parameter space of their models in real-time in an intuitive graphical interface.

You can experiment with miaSimShiny in an [interactive online application](https://gaoyu.shinyapps.io/shiny_rep/).

### Contributions and acknowledgments

You can find us online from [Gitter](https://gitter.im/microbiome/miaverse).

Contributions are very welcome through issues and pull requests at the
[development site](https://github.com/microbiome/miaSim). We a git
flow kind of approach. Development version should be done against the
`main` branch and then merged to `release` for release.
(https://guides.github.com/introduction/flow/)

We are grateful to all
[contributors](https://github.com/microbiome/miaSim/graphs/contributors),
and Emma Gheysen and Karoline Faust for developing the
[microsimR](https://github.com/gheysenemma/microsimR) package that the
miaSim package converts to the SummarizedExperiment universe. This
project is part of [miaverse](microbiome.github.io).

**Kindly cite this work**. For citation details, see R command
  `citation("miaSim")`.


# Code of conduct

Please note that the project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.



[github-watch-badge]: https://img.shields.io/github/watchers/microbiome/miaSim.svg?style=social
[github-watch]: https://github.com/microbiome/miaSim/watchers
[github-star-badge]: https://img.shields.io/github/stars/microbiome/miaSim.svg?style=social
[github-star]: https://github.com/microbiome/miaSim/stargazers
