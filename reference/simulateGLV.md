# Generalized Lotka-Volterra (gLV) simulation

Simulates time series with the generalized Lotka-Volterra model.

## Usage

``` r
simulateGLV(
  n_species,
  names_species = NULL,
  A = NULL,
  x0 = NULL,
  growth_rates = NULL,
  sigma_drift = 0.001,
  sigma_epoch = 0.1,
  sigma_external = 0.3,
  sigma_migration = 0.01,
  epoch_p = 0.001,
  t_external_events = NULL,
  t_external_durations = NULL,
  stochastic = TRUE,
  migration_p = 0.01,
  metacommunity_probability = NULL,
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

- A:

  `Matrix`. Interaction matrix defining the positive and negative
  interactions between n_species. If `NULL`, `randomA(n_species)` is
  used. (Default: `NULL`)

- x0:

  `Numeric scalar`. Indicates the initial abundances of simulated
  species. If `NULL`, `runif(n = n_species, min = 0, max = 1)` is used.
  (Default: `NULL`)

- growth_rates:

  `Numeric scalar`. Indicates the growth rates of simulated species. If
  `NULL`, `runif(n = n_species, min = 0, max = 1)` is used. (Default:
  `NULL`)

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

- stochastic:

  Logical: whether to introduce noise in the simulation. If False,
  sigma_drift, sigma_epoch, and sigma_external are ignored. (default:
  `stochastic = FALSE`)

- migration_p:

  Numeric: the probability/frequency of migration from a metacommunity.
  (default: `migration_p = 0.01`)

- metacommunity_probability:

  Numeric: Normalized probability distribution of the likelihood that
  species from the metacommunity can enter the community during the
  simulation. If NULL, `rdirichlet(1, alpha = rep(1,n_species))` is
  used. (default: `metacommunity_probability = NULL`)

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

`simulateGLV` returns a TreeSummarizedExperiment class object

## Details

Simulates a community time series using the generalized Lotka-Volterra
model, defined as dx/dt = x(b+Ax), where x is the vector of species
abundances, diag(x) is a diagonal matrix with the diagonal values set to
x. A is the interaction matrix and b is the vector of growth rates.

## Examples

``` r

# generate a random interaction matrix
ExampleA <- randomA(n_species = 4, diagonal = -1)

# run the model with default values (only stochastic migration considered)
tse <- simulateGLV(n_species = 4, A = ExampleA)

# run the model with two external disturbances at time points 240 and 480
# with durations equal to 1 (10 time steps when t_step by default is 0.1).
ExampleGLV <- simulateGLV(
    n_species = 4, A = ExampleA,
    t_external_events = c(0, 240, 480), t_external_durations = c(0, 1, 1)
)

# run the model with no perturbation nor migration
set.seed(42)
tse1 <- simulateGLV(
    n_species = 4, A = ExampleA, stochastic = FALSE,
    sigma_migration = 0
)

# run the model with no perturbation nor migration but with measurement error
set.seed(42)
tse2 <- simulateGLV(
    n_species = 4, A = ExampleA, stochastic = FALSE,
    error_variance = 0.001, sigma_migration = 0
)
```
