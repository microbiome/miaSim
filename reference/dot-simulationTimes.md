# Generate simulation times and the indices of time points to return in simulation functions.

Generate simulation times and the indices of time points to return in
simulation functions.

## Usage

``` r
.simulationTimes(t_start = 0, t_end = 1000, t_step = 0.1, t_store = 1000)
```

## Arguments

- t_start:

  `Numeric scalar`. Indicates the initial time of the simulation.
  (Default: `0`)

- t_end:

  `Numeric scalar`. Indicates the final time of the simulation (Default:
  `1000`)

- t_step:

  `Numeric scalar`. Indicates the interval between simulation steps
  (Default: `0.1`)

- t_store:

  `Integer scalar`. Indicates the number of evenly distributed time
  points to keep (Default: `100`)

## Value

lists containing simulation times (t_sys) and the indices to keep.

## Examples

``` r
Time <- .simulationTimes(
    t_start = 0, t_end = 100, t_step = 0.5,
    t_store = 100
)
DefaultTime <- .simulationTimes(t_end = 1000)
```
