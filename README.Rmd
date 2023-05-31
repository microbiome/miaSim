# Microbiome data simulation with miaSim

<!-- badges: start -->

[![Codecov test coverage](https://codecov.io/gh/microbiome/miaSim/branch/master/graph/badge.svg)](https://codecov.io/gh/microbiome/miaSim?branch=master)
[![Gitter](https://badges.gitter.im/microbiome/mia.svg)](https://gitter.im/microbiome/miaverse)
[![Watch on GitHub][github-watch-badge]][github-watch]
[![Star on GitHub][github-star-badge]][github-star]

[![R-CMD-check](https://github.com/microbiome/miaSim/workflows/R-CMD-check/badge.svg)](https://github.com/microbiome/miaSim/actions)
<!-- badges: end -->




## miaSim

This miaSim R/Bioconductor package can be used to simulate
(longitudinal) data for the benchmarking and analysis of quantitative
models of microbial communities.

For installation and use, see the [Getting
started](https://microbiome.github.io/miaSim/articles/miaSim.html)
page.

The package [homepage](https://microbiome.github.io/miaSim/) provides
further tutorials and references for the [implemented
models](https://microbiome.github.io/miaSim/reference/index.html),
which include Hubbell's neutral model, the generalized Lotka-Volterra
model, and the self-organised instability model. The package is based
on the `(Tree)SummarizedExperiment` data container.



### miaSimShiny

The accompanying
[miaSimShiny](https://github.com/gaoyu19920914/miaSimShiny) package
allows users to explore the parameter space of their models in
real-time in an intuitive graphical interface.

You can experiment with miaSimShiny
[online](https://gaoyu.shinyapps.io/shiny_rep/).



### Contributions and acknowledgments

You can find us online from [Gitter](https://gitter.im/microbiome/miaverse).

Contributions are very welcome through issues and pull requests at the
[development site](https://github.com/microbiome/miaSim). We follow a git
flow kind of approach. Development version should be done against the
`main` branch and then merged to `release` for release.
(https://guides.github.com/introduction/flow/)

We are grateful to all
[contributors](https://github.com/microbiome/miaSim/graphs/contributors).


### Citing the package

**Kindly cite this work** as follows:

Gao _et al._ (2023). Methods in Ecology and Evolution. DOI:
[10.1111/2041-210X.14129](https://doi.org/10.1111/2041-210X.14129)

For citation details, see R command `citation("miaSim")`.


# Code of conduct

Please note that the project is released with a [Bioconductor Code of
Conduct](https://bioconductor.github.io/bioc_coc_multilingual/).
By contributing to this project, you agree to abide by its terms.



[github-watch-badge]: https://img.shields.io/github/watchers/microbiome/miaSim.svg?style=social
[github-watch]: https://github.com/microbiome/miaSim/watchers
[github-star-badge]: https://img.shields.io/github/stars/microbiome/miaSim.svg?style=social
[github-star]: https://github.com/microbiome/miaSim/stargazers
