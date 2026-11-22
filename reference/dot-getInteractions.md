# Generate interactions according to five types of interactions and their weights

Generate interactions according to five types of interactions and their
weights

## Usage

``` r
.getInteractions(n_species, weights, connectance)
```

## Arguments

- n_species:

  Integer: defining the dimension of matrix of interaction

- weights:

  Numeric: defining the weights of mutualism, commensalism, parasitism,
  amensalism, and competition in all interspecies interactions.

- connectance:

  Numeric: defining the density of the interaction network. Ranging from
  0 to 1

## Value

A matrix of interactions with all interactions changed according to the
weights and connectance.
