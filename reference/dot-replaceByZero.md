# Replace one element with zero in a list

If the list contains m elements, then lengths of each element must be m,
too. This function is intended to generate a list of x0 (the initial
community) with one missing species, to prepare the parameter
simulations_compare in `estimateAFromSimulations`.

## Usage

``` r
.replaceByZero(input_list)
```

## Arguments

- input_list:

  A list containing m elements, and lengths of each element must be m,
  too.

## Value

A list of same dimension as input_list, but with 0 at specific positions
in the elements of the list.
