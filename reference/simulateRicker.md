# Generate time series with the Ricker model

The Ricker model is a discrete version of the generalized Lotka-Volterra
model and is implemented here as proposed by Fisher and Mehta in PLoS
ONE 2014.

## Usage

``` r
simulateRicker(
  n_species,
  A,
  names_species = NULL,
  x0 = runif(n_species),
  carrying_capacities = runif(n_species),
  error_variance = 0.05,
  explosion_bound = 10^8,
  t_end = 1000,
  norm = FALSE,
  ...
)
```

## Arguments

- n_species:

  Integer: number of species

- A:

  interaction matrix

- names_species:

  Character: names of species. If NULL,
  `paste0("sp", seq_len(n_species))` is used. (default:
  `names_species = NULL`)

- x0:

  `Numeric scalar`. Indicates the initial abundances of simulated
  species. If `NULL`, `runif(n = n_species, min = 0, max = 1)` is used.

- carrying_capacities:

  `Numeric scalar`. Indicates carrying capacities. If NULL,
  `runif(n = n_species, min = 0, max = 1)` is used.

- error_variance:

  `Numeric scalar`. Specifies the variance of measurement error. By
  default it equals to 0, indicating that the result won't contain any
  measurement error. This value should be non-negative. (Default:
  `0.05`)

- explosion_bound:

  `Numeric scalar`. Specifies the boundary for explosion. (Default:
  `10^8`)

- t_end:

  `Integer scalar`. Indicates simulations to be simulated

- norm:

  `Logical scalar`. Whether normalised abundances (proportions in each
  generation) is returned. (Default: `FALSE`)

- ...:

  additional parameters, see
  [`utils`](https://rdrr.io/r/utils/utils-package.html) to know more.

## Value

`simulateRicker` returns a TreeSummarizedExperiment class object

## References

Fisher & Mehta (2014). Identifying Keystone Species in the Human Gut
Microbiome from Metagenomic Timeseries using Sparse Linear Regression.
PLoS One 9:e102451

## Examples

``` r
A <- powerlawA(10, alpha = 1.01)
tse <- simulateRicker(n_species = 10, A, t_end = 100)
```
