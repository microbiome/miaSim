# Stochastic Logistic simulation

Simulates time series with the (stochastic) logistic model

## Usage

``` r
simulateStochasticLogistic(
  n_species,
  names_species = NULL,
  growth_rates = NULL,
  carrying_capacities = NULL,
  death_rates = NULL,
  x0 = NULL,
  sigma_drift = 0.001,
  sigma_epoch = 0.1,
  sigma_external = 0.3,
  sigma_migration = 0.01,
  epoch_p = 0.001,
  t_external_events = NULL,
  t_external_durations = NULL,
  migration_p = 0.01,
  metacommunity_probability = NULL,
  stochastic = TRUE,
  error_variance = 0,
  norm = FALSE,
  t_end = 1000,
  ...
)
```

## Arguments

- n_species:

  Integer: number of species

- names_species:

  Character: names of species. If NULL,
  `paste0("sp", seq_len(n_species))` is used. (default:
  `names_species = NULL`)

- growth_rates:

  `Numeric scalar`. Specifies the growth rates of simulated species. If
  NULL, `runif(n = n_species, min = 0.1, max = 0.2)` is used. (Default:
  `NULL`)

- carrying_capacities:

  `Numeric scalar`. Indicates the max population of species supported in
  the community. If NULL, `runif(n = n_species, min = 1000, max = 2000)`
  is used. (Default: `NULL`)

- death_rates:

  `Numeric scalar`. Indicates the death rates of each species. If NULL,
  `runif(n = n_species, min = 0.0005, max = 0.0025)` is used. (Default:
  `NULL`)

- x0:

  `Numeric scalar`. Indicates the initial abundances of simulated
  species. If NULL, `runif(n = n_species, min = 0.1, max = 10)` is used.
  (Default: `NULL`)

- sigma_drift:

  Numeric: standard deviation of a normally distributed noise applied in
  each time step (t_step) (default: `sigma_drift = 0.001`)

- sigma_epoch:

  Numeric: standard deviation of a normally distributed noise applied to
  random periods of the community composition with frequency defined by
  the epoch_p parameter (default: `sigma_epoch = 0.1`)

- sigma_external:

  Numeric: standard deviation of a normally distributed noise applied to
  user-defined external events/disturbances (default:
  `sigma_external = 0.3`)

- sigma_migration:

  Numeric: standard deviation of a normally distributed variable that
  defines the intensity of migration at each time step (t_step)
  (default: `sigma_migration = 0.01`)

- epoch_p:

  Numeric: the probability/frequency of random periodic changes
  introduced to the community composition (default: `epoch_p = 0.001`)

- t_external_events:

  Numeric: the starting time points of defined external events that
  introduce random changes to the community composition (default:
  `t_external_events = NULL`)

- t_external_durations:

  Numeric: respective duration of the external events that are defined
  in the 't_external_events' (times) and sigma_external (std). (default:
  `t_external_durations = NULL`)

- migration_p:

  Numeric: the probability/frequency of migration from a metacommunity.
  (default: `migration_p = 0.01`)

- metacommunity_probability:

  Numeric: Normalized probability distribution of the likelihood that
  species from the metacommunity can enter the community during the
  simulation. If NULL, `rdirichlet(1, alpha = rep(1,n_species))` is
  used. (default: `metacommunity_probability = NULL`)

- stochastic:

  `Logical scalar`. Whether to introduce noise in the simulation. If
  False, sigma_drift, sigma_epoch, and sigma_external are ignored.
  (Default: `TRUE`)

- error_variance:

  Numeric: the variance of measurement error. By default it equals to 0,
  indicating that the result won't contain any measurement error. This
  value should be non-negative. (default: `error_variance = 0`)

- norm:

  Logical: whether the time series should be returned with the
  abundances as proportions (`norm = TRUE`) or the raw counts (default:
  `norm = FALSE`) (default: `norm = FALSE`)

- t_end:

  Numeric: the end time of the simulationTimes, defining the modeled
  time length of the community. (default: `t_end = 1000`)

- ...:

  additional parameters, see
  [`utils`](https://rdrr.io/r/utils/utils-package.html) to know more.

## Value

`simulateStochasticLogistic` returns a TreeSummarizedExperiment class
object

## Details

The change rate of the species was defined as
`dx/dt = b*x*(1-(x/k))*rN - dr*x`, where b is the vector of growth
rates, x is the vector of initial species abundances, k is the vector of
maximum carrying capacities, rN is a random number ranged from 0 to 1
which changes in each time step, dr is the vector of constant death
rates. Also, the vectors of initial dead species abundances can be set.
The number of species will be set to 0 if the dead species abundances
surpass the alive species abundances.

## Examples

``` r
# Example of logistic model without stochasticity, death rates, or external
# disturbances
set.seed(42)
tse <- simulateStochasticLogistic(
    n_species = 5,
    stochastic = FALSE, death_rates = rep(0, 5)
)

# Adding a death rate
set.seed(42)
tse1 <- simulateStochasticLogistic(
    n_species = 5,
    stochastic = FALSE, death_rates = rep(0.01, 5)
)

# Example of stochastic logistic model with measurement error
set.seed(42)
tse2 <- simulateStochasticLogistic(
    n_species = 5,
    error_variance = 1000
)

# example with all the initial parameters defined by the user
set.seed(42)
tse3 <- simulateStochasticLogistic(
    n_species = 2,
    names_species = c("species1", "species2"),
    growth_rates = c(0.2, 0.1),
    carrying_capacities = c(1000, 2000),
    death_rates = c(0.001, 0.0015),
    x0 = c(3, 0.1),
    sigma_drift = 0.001,
    sigma_epoch = 0.3,
    sigma_external = 0.5,
    sigma_migration = 0.002,
    epoch_p = 0.001,
    t_external_events = c(100, 200, 300),
    t_external_durations = c(0.1, 0.2, 0.3),
    migration_p = 0.01,
    metacommunity_probability = miaSim::rdirichlet(1, alpha = rep(1, 2)),
    stochastic = TRUE,
    error_variance = 0,
    norm = FALSE, # TRUE,
    t_end = 400,
    t_start = 0, t_step = 0.01,
    t_store = 1500
)
```
