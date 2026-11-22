# Hubbell's neutral model simulation

Neutral species abundances simulation according to the Hubbell model.

## Usage

``` r
simulateHubbell(
  n_species,
  M,
  carrying_capacity = 1000,
  k_events = 10,
  migration_p = 0.02,
  t_skip = 0,
  t_end,
  norm = FALSE
)
```

## Arguments

- n_species:

  `Integer scalar`. Specifies the amount of different species initially
  in the local community.

- M:

  `Integer scalar`. Specifies the amount of different species in the
  metacommunity, including those of the local community

- carrying_capacity:

  `Integer scalar`. Indicates the fixed amount of individuals in the
  local community. (Default: `1000`)

- k_events:

  `Inteer scalar`. Indicates the fixed amount of deaths of local
  community individuals in each generation (Default: `10`)

- migration_p:

  `Numeric scalar`. The immigration rate; specifies the probability that
  a death in the local community is replaced by a migrant of the
  metacommunity rather than by the birth of a local community member
  (Default: `0.02`)

- t_skip:

  `Integer scalar`. Indicates the number of generations that should not
  be included in the outputted species abundance matrix. (Default: `0`)

- t_end:

  `Integer scalar`. Indicates the number of simulations to be simulated

- norm:

  `Logical scalar`. Whether the time series should be returned with the
  abundances as proportions (`norm = TRUE`) or the raw counts. (Default:
  `FALSE`)

## Value

`simulateHubbell` returns a TreeSummarizedExperiment class object

## References

Rosindell, James et al. "The unified neutral theory of biodiversity and
biogeography at age ten." Trends in ecology & evolution vol. 26,7
(2011).

## Examples

``` r
tse <- simulateHubbell(
    n_species = 8, M = 10, carrying_capacity = 1000, k_events = 50,
    migration_p = 0.02, t_end = 100
)
```
