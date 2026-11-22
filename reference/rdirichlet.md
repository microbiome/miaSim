# Generate dirichlet random deviates

Generate dirichlet random deviates

## Usage

``` r
rdirichlet(n, alpha)
```

## Arguments

- n:

  Number of random vectors to generate.

- alpha:

  Vector containing shape parameters.

## Value

a vector containing the Dirichlet density

## Examples

``` r
dirichletExample <- rdirichlet(1, c(1, 2, 3))
```
