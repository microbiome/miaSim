# Consumer-resource model simulation

Simulates time series with the consumer-resource model.

## Usage

``` r
simulateConsumerResource(
  n_species,
  n_resources,
  names_species = NULL,
  names_resources = NULL,
  E = NULL,
  x0 = NULL,
  resources = NULL,
  resources_dilution = NULL,
  growth_rates = NULL,
  monod_constant = NULL,
  sigma_drift = 0.001,
  sigma_epoch = 0.1,
  sigma_external = 0.3,
  sigma_migration = 0.01,
  epoch_p = 0.001,
  t_external_events = NULL,
  t_external_durations = NULL,
  stochastic = FALSE,
  migration_p = 0.01,
  metacommunity_probability = NULL,
  error_variance = 0,
  norm = FALSE,
  t_end = 1000,
  trophic_priority = NULL,
  inflow_rate = 0,
  outflow_rate = 0,
  volume = 1000,
  ...
)
```

## Arguments

- n_species:

  Integer: number of species

- n_resources:

  Integer: number of resources

- names_species:

  Character: names of species. If NULL,
  `paste0("sp", seq_len(n_species))` is used. (default:
  `names_species = NULL`)

- names_resources:

  Character: names of resources. If NULL,
  `paste0("res", seq_len(n_resources))` is used.

- E:

  `Matrix`. Defines the efficiency of resource-to-biomass conversion
  (positive values) and the relative conversion of metabolic by-products
  (negative values). If `NULL`, `randomE(n_species, n_resources)` is
  used. (Default: `NULL`)

- x0:

  `Numeric scalar`. Specifies the initial abundances of simulated
  species. If `NULL`, `runif(n = n_species, min = 0.1, max = 10)` is
  used. (Default: `NULL`)

- resources:

  `Numeric scalar`. Specifies the initial concentrations of resources.
  If `NULL`, `runif(n = n_resources, min = 1, max = 100)` is used.
  (Default: `NULL`)

- resources_dilution:

  `Numeric scalar`. Specifies the concentrations of resources in the
  continuous inflow (applicable when inflow_rate \> 0). If `NULL`,
  `resources` is used. (Default: `NULL`)

- growth_rates:

  `Numeric vector`. Specifies the maximum growth rates(mu) of species.
  If `NULL`, `rep(1, n_species)` is used. (Default: `NULL`)

- monod_constant:

  `Matrix`. Specifies the constant of additive monod growth of n_species
  consuming n_resources. If `NULL`,
  `matrix(rgamma(n = n_species*n_resources, shape = 50*max(resources), rate = 1), nrow = n_species)`
  is used. (Default: `NULL`)

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

- trophic_priority:

  `Matrix`. Defines the orders of resources to be consumed by each
  species. If `NULL`, by default, this feature won't be turned on, and
  species will consume all resources simultaneously to grow. The
  dimension should be identical to matrix E. (Default: `NULL`)

- inflow_rate, :

  outflow_rate `Numeric scalar`. The inflow of a culture process. By
  default, inflow_rate and is 0, indicating a batch culture process.
  When larger than 0, we can simulate a continuous culture(e.g.
  chemostat).

- outflow_rate:

  `Numeric scalar`. outflow rate of a culture process By default,
  outflow_rate is 0, indicating a batch culture process. When larger
  than 0, we can simulate a continuous culture(e.g. chemostat).

- volume:

  `Numeric scalar`. Indicates the volume of the continuous cultivation.
  This parameter is important for simulations where inflow_rate or
  outflow_rate are not 0. (Default: `1000`)

- ...:

  additional parameters, see
  [`utils`](https://rdrr.io/r/utils/utils-package.html) to know more.

## Value

an TreeSummarizedExperiment class object

## Examples

``` r
n_species <- 2
n_resources <- 4
tse <- simulateConsumerResource(
    n_species = n_species,
    n_resources = n_resources
)

if (FALSE) { # \dontrun{
# example with user-defined values (names_species, names_resources, E, x0,
# resources, growth_rates, error_variance, t_end, t_step)

ExampleE <- randomE(
    n_species = n_species, n_resources = n_resources,
    mean_consumption = 3, mean_production = 1, maintenance = 0.4
)
ExampleResources <- rep(100, n_resources)
tse1 <- simulateConsumerResource(
    n_species = n_species,
    n_resources = n_resources, names_species = letters[seq_len(n_species)],
    names_resources = paste0("res", LETTERS[seq_len(n_resources)]), E = ExampleE,
    x0 = rep(0.001, n_species), resources = ExampleResources,
    growth_rates = runif(n_species),
    error_variance = 0.01,
    t_end = 5000,
    t_step = 1
)

# example with trophic levels
n_species <- 10
n_resources <- 15
ExampleEfficiencyMatrix <- randomE(
    n_species = 10, n_resources = 15,
    trophic_levels = c(6, 3, 1),
    trophic_preferences = list(
        c(rep(1, 5), rep(-1, 5), rep(0, 5)),
        c(rep(0, 5), rep(1, 5), rep(-1, 5)),
        c(rep(0, 10), rep(1, 5))
    )
)

ExampleResources <- c(rep(500, 5), rep(200, 5), rep(50, 5))
tse2 <- simulateConsumerResource(
    n_species = n_species,
    n_resources = n_resources,
    names_species = letters[1:n_species],
    names_resources = paste0(
        "res", LETTERS[1:n_resources]
    ),
    E = ExampleEfficiencyMatrix,
    x0 = rep(0.001, n_species),
    resources = ExampleResources,
    growth_rates = rep(1, n_species),
    # error_variance = 0.001,
    t_end = 5000, t_step = 1
)

# example with trophic priority
n_species <- 4
n_resources <- 6
ExampleE <- randomE(
    n_species = n_species,
    n_resources = n_resources,
    mean_consumption = n_resources,
    mean_production = 0
)
ExampleTrophicPriority <- t(apply(
    matrix(sample(n_species * n_resources),
        nrow = n_species
    ),
    1, order
))
# make sure that for non-consumables resources for each species,
# the priority is 0 (smaller than any given priority)
ExampleTrophicPriority <- (ExampleE > 0) * ExampleTrophicPriority
tse3 <- simulateConsumerResource(
    n_species = n_species,
    n_resources = n_resources,
    E = ExampleE,
    trophic_priority = ExampleTrophicPriority,
    t_end = 2000
)
} # }
```
