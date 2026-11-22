# Get the interspecies interaction matrix A using leave-one-out method

generate matrix A from the comparisons between simulations with one
absent species and a simulation with complete species (leave-one-out)

## Usage

``` r
.estimateAFromSimulations(
  simulations,
  simulations2,
  n_instances = 1,
  t_end = NULL,
  scale_off_diagonal = 0.1,
  diagonal = -0.5,
  connectance = 0.2
)
```

## Arguments

- simulations:

  A list of simulation(s) with complete species

- simulations2:

  A list of simulation(s), each with one absent species

- n_instances:

  Integer: number of instances to generate (default: `n_instances = 1`)

- t_end:

  Numeric: end time of the simulation. If not identical with t_end in
  params_list, then it will overwrite t_end in each simulation (default:
  `t_end = 1000`)

- scale_off_diagonal:

  Numeric: scale of the off-diagonal elements compared to the diagonal.
  Same to the parameter in function `randomA`. (default:
  `scale_off_diagonal = 0.1`)

- diagonal:

  Values defining the strength of self-interactions. Input can be a
  number (will be applied to all species) or a vector of length
  n_species. Positive self-interaction values lead to exponential
  growth. Same to the parameter in function `randomA`. (default:
  `diagonal = -0.5`)

- connectance:

  Numeric frequency of inter-species interactions. i.e. proportion of
  non-zero off-diagonal terms. Should be in the interval 0 \<=
  connectance \<= 1. Same to the parameter in function `randomA`.
  (default: `connectance = 0.2`)

## Value

a matrix A with dimensions (n_species x n_species) where n_species
equals to the number of elements in simulations2
