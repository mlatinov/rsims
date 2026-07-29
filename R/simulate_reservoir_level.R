#' Simulate Reservoir Water Levels with a Bayesian Structural Time Series
#'
#' @description
#' Simulates daily reservoir water levels using a simple Bayesian Structural
#' Time Series (BSTS) model composed of:
#'
#' * a local-level state that evolves as a random walk,
#' * a piecewise-constant regime component representing seasonal hydrological
#'   periods (for example wet and dry seasons),
#' * Gaussian observation noise.
#'
#' The latent state evolves according to
#'
#' \deqn{
#' \mu_t = \mu_{t-1} + \epsilon_t,
#' \qquad
#' \epsilon_t \sim Normal(0,\sigma_{state}).
#' }
#'
#' The observed reservoir level is
#'
#' \deqn{
#' y_t = \mu_t + s_t + \varepsilon_t,
#' }
#'
#' where
#'
#' * \eqn{\mu_t} is the latent reservoir level,
#' * \eqn{s_t} is the seasonal/regime effect,
#' * \eqn{\varepsilon_t \sim Normal(0,\sigma_{obs})} is observation noise.
#'
#' This simulation is useful for demonstrating local-level state-space models,
#' Bayesian Structural Time Series (BSTS), Kalman filtering, and forecasting.
#'
#' @param T Integer.
#' Number of time points to simulate.
#'
#' @param inital_reservoir_levels Numeric.
#' Initial latent reservoir level at the beginning of the series.
#'
#' @param inital_reservoir_levels_sd Numeric.
#' Standard deviation of the latent state innovations controlling how much
#' the underlying reservoir level changes from one time point to the next.
#'
#' @param reservoir_regime_lenght Integer.
#' Number of consecutive observations belonging to each seasonal regime.
#'
#' @param regime_effects Numeric vector.
#' Additive effect associated with each regime. Positive values increase the
#' expected reservoir level whereas negative values decrease it.
#'
#' @param sigma_reservoir_levels Numeric.
#' Standard deviation of the Gaussian observation error.
#'
#' @return
#' A data frame containing:
#'
#' * `reservoir_level` — observed reservoir level.
#' * `day` — time index.
#'
#' @details
#' Internally this function combines three simulation helpers:
#'
#' * `simulate_local_level()` generates the latent random-walk state.
#' * `simulate_regime()` creates piecewise-constant seasonal regimes.
#' * `simulate_bsts()` combines the state, seasonal component and observation
#'   model into the final BSTS simulation.
#'
#' @examples
#' reservoir <- simulate_reservoir_level()
#'
#' plot(
#'   reservoir$day,
#'   reservoir$reservoir_level,
#'   type = "l",
#'   xlab = "Day",
#'   ylab = "Reservoir level"
#' )
#'
#' @export
simulate_reservoir_level <- function(
  T = 200,
  initial_reservoir_levels    = 50,
  initial_reservoir_levels_sd = 0.5,
  reservoir_regime_length    = 50,
  regime_effects             = c(5, 4, 2, -2 ),
  sigma_reservoir_levels     = 0.5

){

  # Simulate the reservoir baseline level for the period T
  level  <- simulate_local_level(T, initial = initial_reservoir_levels, sigma = initial_reservoir_levels_sd)

  # Simulate Sesonalal regimes every J months
  regime <- simulate_regime(
    T = T,
    block_length = reservoir_regime_length,
    effects      = regime_effects
  )

  # Simulate BSST type model that follows the Reservoir Levels
  reservoir_levels <- simulate_bsts(
    T     = T,
    state = level$state,
    seasonal = regime,
    family = "gaussian",
    sigma_obs = sigma_reservoir_levels
  )

  # Combine and return a dataframe with the Resevoir Levels and days 
  sim_data <- data.frame(
    reservoir_level = reservoir_levels$y,
    day             = 1:T
  )
  return(sim_data)
}