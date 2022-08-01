#' @param sigma_drift Numeric: standard deviation of a normally distributed
#' noise applied in each time step (t_step)
#' (default: \code{sigma_drift = 0.001})
#' @param sigma_epoch Numeric: standard deviation of a normally distributed
#' noise applied to random periods of the community composition with frequency
#' defined by the epoch_p parameter
#' (default: \code{sigma_epoch = 0.1})
#' @param sigma_external Numeric: standard deviation of a normally distributed
#' noise applied to user-defined external events/disturbances
#' (default: \code{sigma_external = 0.3})
#' @param sigma_migration Numeric: standard deviation of a normally distributed
#' variable that defines the intensity of migration at each time step (t_step)
#' (default: \code{sigma_migration = 0.01})
#' @param epoch_p Numeric: the probability/frequency of random periodic
#' changes introduced to the community composition
#' (default: \code{epoch_p = 0.001})
#' @param t_external_events Numeric: the starting time points of defined
#' external events that introduce random changes to the community composition
#' (default: \code{t_external_events = NULL})
#' @param t_external_durations Numeric: respective duration of the external
#' events that are defined in the 't_external_events' (times) and
#' sigma_external (std).
#' (default: \code{t_external_durations = NULL})
#' @param stochastic Logical: whether to introduce noise in the simulation.
#' If False, sigma_drift, sigma_epoch, and sigma_external are ignored.
#' (default: \code{stochastic = FALSE})
