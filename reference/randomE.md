# Generate random efficiency matrix

Generate random efficiency matrix for consumer resource model from
Dirichlet distribution, where positive efficiencies indicate the
consumption of resources, whilst negatives indicate that the species
would produce the resource.

## Usage

``` r
randomE(
  n_species,
  n_resources,
  names_species = NULL,
  names_resources = NULL,
  mean_consumption = n_resources/4,
  mean_production = n_resources/6,
  maintenance = 0.5,
  trophic_levels = NULL,
  trophic_preferences = NULL,
  exact = FALSE
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

- mean_consumption:

  `Numeric scalar`. Specifies the mean number of resources consumed by
  each species drawn from a poisson distribution (Default:
  `n_resources/4`)

- mean_production:

  `Numeric scalar`. Specifies the mean number of resources produced by
  each species drawn from a poisson distribution (Default:
  `n_resources/6`)

- maintenance:

  `Numeric scalar`. Specifies the proportion of resources that cannot be
  converted into products between 0~1 the proportion of resources used
  to maintain the living of microorganisms. 0 means all the resources
  will be used for the reproduction of microorganisms, and 1 means all
  the resources would be used to maintain the living of organisms and no
  resources would be left for their growth(reproduction). (Default:
  `0.5`)

- trophic_levels:

  `Integer scalar`. Indicates the number of species in microbial trophic
  levels. If NULL, by default, microbial trophic levels would not be
  considered. (Default: `NULL`)

- trophic_preferences:

  `List`. Indicates the preferred resources and productions of each
  trophic level. Positive values indicate the consumption of resources,
  whilst negatives indicate that the species would produce the resource.
  (Default: `NULL`)

- exact:

  `Logical scalar`. Whether to set the number of consumption/production
  to be exact as mean_consumption/mean_production or to set them using a
  Poisson distribution. (Default: `FALSE`) If
  `length(trophic_preferences)` is smaller than
  `length(trophic_levels)`, then NULL values would be appended to lower
  trophic levels. If NULL, by default, the consumption preference will
  be defined randomly. (Default: `trophic_preferences = NULL`)

## Value

`randomE` returns a matrix E with dimensions (n_species x n_resources),
and each row represents a species.

## Examples

``` r
# example with minimum parameters
ExampleEfficiencyMatrix <- randomE(n_species = 5, n_resources = 12)

# examples with specific parameters
ExampleEfficiencyMatrix <- randomE(
    n_species = 3, n_resources = 6,
    names_species = letters[1:3],
    names_resources = paste0("res", LETTERS[1:6]),
    mean_consumption = 3, mean_production = 1
)
ExampleEfficiencyMatrix <- randomE(
    n_species = 3, n_resources = 6,
    maintenance = 0.4
)
ExampleEfficiencyMatrix <- randomE(
    n_species = 3, n_resources = 6,
    mean_consumption = 3, mean_production = 1, maintenance = 0.4
)

# examples with microbial trophic levels
ExampleEfficiencyMatrix <- randomE(
    n_species = 10, n_resources = 15,
    trophic_levels = c(6, 3, 1),
    trophic_preferences = list(
        c(rep(1, 5), rep(-1, 5), rep(0, 5)),
        c(rep(0, 5), rep(1, 5), rep(-1, 5)),
        c(rep(0, 10), rep(1, 5))
    )
)
ExampleEfficiencyMatrix <- randomE(
    n_species = 10, n_resources = 15,
    trophic_levels = c(6, 3, 1),
    trophic_preferences = list(c(rep(1, 5), rep(-1, 5), rep(0, 5)), NULL, NULL)
)
ExampleEfficiencyMatrix <- randomE(
    n_species = 10, n_resources = 15,
    trophic_levels = c(6, 3, 1)
)
```
