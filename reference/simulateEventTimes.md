# Generate a vector of event times

Generate a vector of event times

## Usage

``` r
simulateEventTimes(t_events = NULL, t_duration = NULL, t_end = 1000, ...)
```

## Arguments

- t_events:

  Numeric vector; starting time of the events

- t_duration:

  Numeric vector; duration of the events

- t_end:

  Numeric: end time of the simulation

- ...:

  : additional parameters to pass to simulationTimes, including t_start,
  t_step, and t_store.

## Value

A vector of time points in the simulation

## Examples

``` r
tEvent <- simulateEventTimes(
    t_events = c(10, 50, 100),
    t_duration = c(1, 2, 3),
    t_end = 100,
    t_store = 100,
    t_step = 1
)
```
