# Generate pairs of interactions according to interaction types

A helper function to be used in combination with .getInteractions()

## Usage

``` r
.applyInterType(I, pair, interType)
```

## Arguments

- I:

  Matrix: defining the interaction between each pair of species

- pair:

  Numeric: a vector with a length of 2, indicating the 2 focusing
  species in the process of applying the interaction types

- interType:

  Character: one of 'mutualism', 'commensalism', 'parasitism',
  'amensalism', or 'competition'. Defining the interaction type

## Value

A matrix of interaction types with one pair changed
