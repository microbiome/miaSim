# Self-Organised Instability model (SOI) simulation

Generate time-series with The Self-Organised Instability (SOI) model.
Implements a K-leap method for accelerating stochastic simulation.

## Usage

``` r
simulateSOI(
  n_species,
  x0 = NULL,
  names_species = NULL,
  carrying_capacity = 1000,
  A = NULL,
  k_events = 5,
  t_end = 1000,
  metacommunity_probability = runif(n_species, min = 0.1, max = 0.8),
  death_rates = runif(n_species, min = 0.01, max = 0.08),
  norm = FALSE
)
```

## Arguments

- n_species:

  Integer: number of species

- x0:

  `Numeric scalar`. Specifies initial community abundances If `NULL`,
  based on migration rates. (Default: `NULL`)

- names_species:

  Character: names of species. If NULL,
  `paste0("sp", seq_len(n_species))` is used. (default:
  `names_species = NULL`)

- carrying_capacity:

  `Integer scalar`. Indicates community size, number of available sites
  (individuals). (Default: `1000`)

- A:

  `Matrix`. Defines the positive and negative interactions between
  n_species. If `NULL`, `powerlawA(n_species)` is used. (Default:
  `NULL`)

- k_events:

  `Integer scalar`. Indicates the number of transition events that are
  allowed to take place during one leap. Higher values reduce runtime,
  but also accuracy of the simulation. (Default: `5`).

- t_end:

  `Numeric scalar`. Specifies the end time of the simulation, defining
  the modeled time length of the community. (Default: `1000`)

- metacommunity_probability:

  `Numeric scalar`: Indicates the normalized probability distribution of
  the likelihood that species from the metacommunity can enter the
  community during the simulation. (Default:
  `runif(n_species, min = 0.1, max = 0.8)`)

- death_rates:

  `Numeric scalar`. Indicates the death rates of each species. (Default:
  `runif(n_species, min = 0.01, max = 0.08)`)

- norm:

  `Logical scalar`. Whether the time series should be returned with the
  abundances as proportions (`norm = TRUE`) or the raw counts. (Default:
  `FALSE`)

## Value

`simulateSOI` returns a TreeSummarizedExperiment class object

## Examples

``` r
# Generate interaction matrix
A <- miaSim::powerlawA(10, alpha = 1.2)
# Simulate data from the SOI model
tse <- simulateSOI(
    n_species = 10, carrying_capacity = 1000, A = A,
    k_events = 5, x0 = NULL, t_end = 150, norm = TRUE
)
```
