# miaSim: Microbiome Data Simulation

## Introduction

`miaSim` implements tools for microbiome data simulation based on
varying ecological modeling assumptions. These can be used to simulate
species abundance matrices, including time series. Detailed function
documentation is available at the [function
reference](https://microbiome.github.io/miaSim/reference/index.html)

The miaSim package supports the R/Bioconductor multi-assay framework.
For more information on operating with this data format in microbial
ecology, see the [online tutorial](https://microbiome.github.io/OMA).

### Installation

The stable Bioconductor release version can be installed as follows.

\
`if`` ``(``!`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"BiocManager"``, quietly ``=`` ``TRUE``)``)`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"BiocManager"``)`\
`if`` ``(``!`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"miaSim"``, quietly ``=`` ``TRUE``)``)``    `\
`    ``BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"miaSim"``)`

The experimental Bioconductor devel version can be installed as follows.

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    # The following initializes usage of Bioc devel
    BiocManager::install(version='devel')
    BiocManager::install("miaSim")

Load the library

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`miaSim`](https://github.com/microbiome/miaSim)`)`

### Examples

#### Generate species interaction matrices for the models

Some of the models rely on interaction matrices that represents
interaction heterogeneity between species. The interaction matrix can be
generated with different distributional assumptions.

Generate interactions from normal distribution:

\
`A_normal`` ``<-`` `[`powerlawA`](../reference/powerlawA.md)`(``n_species ``=`` ``4``, alpha ``=`` ``3``)`

Generate interactions from uniform distribution:

\
`A_uniform`` ``<-`` `[`randomA`](../reference/randomA.md)`(``n_species ``=`` ``10``,`\
`                 diagonal ``=`` ``-``0.4``,`\
`                     connectance ``=`` ``0.5``,`\
`             interactions ``=`` `[`runif`](https://rdrr.io/r/stats/Uniform.html)`(``n ``=`` ``10``^``2``, min ``=`` ``-``0.8``, max ``=`` ``0.8``)``)`

#### Hubbell model

Hubbell Neutral simulation model characterizes diversity and relative
abundance of species in ecological communities assuming migration,
births and deaths but no interactions. Losses become replaced by
migration or birth.

\
`tse_hubbell`` ``<-`` `[`simulateHubbell`](../reference/simulateHubbell.md)`(``n_species ``=`` ``8``,`\
`                               M ``=`` ``10``,`\
`                   carrying_capacity ``=`` ``1000``,`\
`                               k_events ``=`` ``50``,`\
`                   migration_p ``=`` ``0.02``,`\
`                   t_end ``=`` ``100``)`

One can also simulate parameters for the Hubbell model.

\
`params_hubbell`` ``<-`` `[`simulateHubbellRates`](../reference/simulateHubbellRates.md)`(``x0 ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``,``5``,``10``)``,`\
`    migration_p ``=`` ``0.1``, metacommunity_probability ``=`` ``NULL``, k_events ``=`` ``1``, `\
`    growth_rates ``=`` ``NULL``, norm ``=`` ``FALSE``, t_end``=``1000``)`

#### Stochastic logistic model

Stochastic logistic model is used to determine dead and alive counts in
community.

\
`tse_logistic`` ``<-`` `[`simulateStochasticLogistic`](../reference/simulateStochasticLogistic.md)`(``n_species ``=`` ``5``)`

#### Self-Organised Instability (SOI)

The Self-Organised Instability (SOI) model generates time series for
communities and accelerates stochastic simulation.

\
`tse_soi`` ``<-`` `[`simulateSOI`](../reference/simulateSOI.md)`(``n_species ``=`` ``4``, carrying_capacity ``=`` ``1000``,`\
`                       A ``=`` ``A_normal``, k_events``=``5``,`\
`               x0 ``=`` ``NULL``,t_end ``=`` ``150``, norm ``=`` ``TRUE``)`

#### Consumer-resource model

The consumer resource model requires the `randomE` function. This
returns a matrix containing the production rates and consumption rates
of each species. The resulting matrix is used as a determination of
resource consumption efficiency.

\
`# Consumer-resource model as a TreeSE object`\
`tse_crm`` ``<-`` `[`simulateConsumerResource`](../reference/simulateConsumerResource.md)`(``n_species ``=`` ``2``,`\
`                                    n_resources ``=`` ``4``,`\
`                                    E ``=`` `[`randomE`](../reference/randomE.md)`(``n_species ``=`` ``2``, n_resources ``=`` ``4``)``)`

You could visualize the simulated dynamics using tools from the
[miaTime](https://microbiome.github.io/miaTime/) package.

#### Generalized Lotka-Volterra (gLV)

The generalized Lotka-Volterra simulation model generates time-series
assuming microbial population dynamics and interaction.

\
`tse_glv`` ``<-`` `[`simulateGLV`](../reference/simulateGLV.md)`(``n_species ``=`` ``4``,`\
`                       A ``=`` ``A_normal``,`\
`               t_start ``=`` ``0``, `\
`                       t_store ``=`` ``1000``,`\
`               stochastic ``=`` ``FALSE``,`\
`               norm ``=`` ``FALSE``)`

#### Ricker model

Ricker model is a discrete version of the gLV:

\
`tse_ricker`` ``<-`` `[`simulateRicker`](../reference/simulateRicker.md)`(``n_species``=``4``, A ``=`` ``A_normal``, t_end``=``100``, norm ``=`` ``FALSE``)`

The number of species specified in the interaction matrix must be the
same as the species used in the models.

### Data containers

The simulated data sets are returned as `TreeSummarizedExperiment`
objects. This provides access to a broad range of tools for microbiome
analysis that support this format (see
[microbiome.github.io](http://microbiome.github.io)). More examples on
the object manipulation and analysis can be found at [OMA Online
Manual](https://microbiome.github.io/OMA).

For instance, to plot population density we can use the `miaViz`
package:

### Case studies

Source code for replicating the published case studies using the miaSim
package ([Gao et al. 2023](https://doi.org/10.1111/2041-210X.14129)) is
available in
[Github](https://github.com/microbiome/miaSim/tree/main/inst/extdata/phyloseq)
(based on the phyloseq data container).

### Related work

- [micodymora](https://github.com/OSS-Lab/micodymora) Python package for
  microbiome simulation

## Session info

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] miaSim_1.15.1                   TreeSummarizedExperiment_2.21.0
    ##  [3] Biostrings_2.81.9               XVector_0.53.0                 
    ##  [5] SingleCellExperiment_1.35.2     SummarizedExperiment_1.43.0    
    ##  [7] Biobase_2.73.2                  GenomicRanges_1.65.4           
    ##  [9] Seqinfo_1.3.2                   IRanges_2.47.5                 
    ## [11] S4Vectors_0.51.10               BiocGenerics_0.59.13           
    ## [13] generics_0.1.4                  MatrixGenerics_1.25.0          
    ## [15] matrixStats_1.5.0               BiocStyle_2.41.0               
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] xfun_0.62           bslib_0.12.0        poweRlaw_1.0.0     
    ##  [4] htmlwidgets_1.6.4   lattice_0.23-1      yulab.utils_0.2.6  
    ##  [7] vctrs_0.7.3         tools_4.6.1         parallel_4.6.1     
    ## [10] tibble_3.3.1        pkgconfig_2.0.3     Matrix_1.7-6       
    ## [13] desc_1.4.3          lifecycle_1.0.5     compiler_4.6.1     
    ## [16] treeio_1.37.1       textshaping_1.0.5   codetools_0.2-20   
    ## [19] htmltools_0.5.9     sass_0.4.10         yaml_2.3.12        
    ## [22] lazyeval_0.2.3      pkgdown_2.2.1       pillar_1.11.1      
    ## [25] crayon_1.5.3        jquerylib_0.1.4     tidyr_1.3.2        
    ## [28] BiocParallel_1.47.0 DelayedArray_0.39.8 cachem_1.1.0       
    ## [31] abind_1.4-8         nlme_3.1-171        tidyselect_1.2.1   
    ## [34] digest_0.6.39       purrr_1.2.2         dplyr_1.2.1        
    ## [37] bookdown_0.48       fastmap_1.2.0       grid_4.6.1         
    ## [40] cli_3.6.6           SparseArray_1.13.5  magrittr_2.0.5     
    ## [43] S4Arrays_1.13.3     ape_5.8-1           rappdirs_0.3.4     
    ## [46] rmarkdown_2.32      otel_0.2.0          deSolve_1.42       
    ## [49] ragg_1.5.2          evaluate_1.0.5      knitr_1.52         
    ## [52] rlang_1.3.0         Rcpp_1.1.2          glue_1.8.1         
    ## [55] tidytree_0.4.9      BiocManager_1.30.27 jsonlite_2.0.0     
    ## [58] R6_2.6.1            systemfonts_1.3.2   fs_2.1.0
