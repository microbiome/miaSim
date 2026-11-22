# Interaction matrix with Power-Law network adjacency matrix

N is the an Interspecific Interaction matrix with values drawn from a
normal distribution H the interaction strength heterogeneity drawn from
a power-law distribution with the parameter alpha, and G the adjacency
matrix of with out-degree that reflects the heterogeneity of the
powerlaw. A scaling factor s may be used to constrain the values of the
interaction matrix to be within a desired range. Diagonal elements of A
are defined by the parameter d.

## Usage

``` r
powerlawA(n_species, alpha = 3, stdev = 1, s = 0.1, d = -1, symmetric = FALSE)
```

## Arguments

- n_species:

  `Integer scalar`. Indicates the number of species.

- alpha:

  `Numeric scalar`. Specifies the power-law distribution. Should be
  \> 1. Larger values will give lower interaction strength
  heterogeneity, whereas values closer to 1 give strong heterogeneity in
  interaction strengths between the species. In other words, values of
  alpha close to 1 will give Strongly Interacting Species (SIS).
  (Default: `3.0`)

- stdev:

  `Numeric scalar`. Specifies the standard deviation of the normal
  distribution with mean 0 from which the elements of the nominal
  interspecific interaction matrix N are drawn. (Default: `1`)

- s:

  `Numeric scalar`. Specifies the scaling with which the final global
  interaction matrix A is multiplied. (Default: `0.1`)

- d:

  `Numeric scalar`. Diagonal values, indicating self-interactions (use
  negative values for stability). (Default: `1.0`)

- symmetric:

  `Logical scalar`. Whether a symmetric interaction matrix is returned.
  (Default: `FALSE`)

## Value

The interaction matrix A with dimensions (n_species x n_species)

## References

Gibson TE, Bashan A, Cao HT, Weiss ST, Liu YY (2016) On the Origins and
Control of Community Types in the Human Microbiome. PLOS Computational
Biology 12(2): e1004688. https://doi.org/10.1371/journal.pcbi.1004688

## Examples

``` r
# Low interaction heterogeneity
A_low <- powerlawA(n_species = 10, alpha = 3)
# Strong interaction heterogeneity
A_strong <- powerlawA(n_species = 10, alpha = 1.01)
```
