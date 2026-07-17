#' Simulate Solar Panel Efficiency Measurements
#'
#' @description
#' ## Scenario
#'
#' A renewable energy company monitors multiple solar installations over
#' several days to understand how environmental conditions affect the efficiency
#' of photovoltaic panels.
#'
#' For each daily measurement you observe:
#'
#' * air temperature,
#' * a dust accumulation index,
#' * the installation where the measurement was collected, and
#' * the observed panel efficiency.
#'
#' The relationship between temperature and efficiency is assumed to be
#' nonlinear. Rather than specifying a parametric equation, the simulator
#' generates the temperature effect using radial basis functions (Gaussian
#' bumps) whose weights define the shape of the curve.
#'
#' Different installations are also allowed to have different baseline
#' efficiencies and different sensitivities to dust accumulation. These two
#' installation-level effects are sampled jointly from a correlated multivariate
#' normal distribution, allowing installations with higher baseline efficiency
#' to systematically differ in how strongly dust affects performance.
#'
#' The challenge is to recover the nonlinear temperature response while
#' accounting for hierarchical variation among installations.
#'
#' ## Data Generating Model
#'
#' For measurement \eqn{i} from installation \eqn{j},
#'
#' \deqn{
#' Efficiency_i \sim Beta(\mu_i\kappa,\ (1-\mu_i)\kappa)
#' }
#'
#' where
#'
#' \deqn{
#' logit(\mu_i)=
#' \alpha_j
#' +\beta_{Dust,j}Dust_i
#' +f(Temperature_i).
#' }
#'
#' The nonlinear temperature effect is represented as
#'
#' \deqn{
#' f(x)=
#' \sum_{k=1}^{K}
#' w_kB_k(x),
#' }
#'
#' where \eqn{B_k(x)} are radial Gaussian basis functions.
#'
#' Installation-specific effects are generated jointly as
#'
#' \deqn{
#' \begin{pmatrix}
#' \alpha_j\\
#' \beta_{Dust,j}
#' \end{pmatrix}
#' \sim
#' MVN
#' \left(
#' \begin{pmatrix}
#' \bar{\alpha}\\
#' \bar{\beta}
#' \end{pmatrix},
#' \Sigma
#' \right).
#' }
#'
#' This simulator therefore combines a nonlinear spline effect with correlated
#' varying intercepts and varying slopes.
#'
#' @param num_installations Integer. Number of solar installations.
#'
#' @param num_days Integer. Number of daily observations recorded per
#' installation.
#'
#' @param mean_temperature Mean daily air temperature (°C).
#'
#' @param sd_temperature Standard deviation of daily temperatures.
#'
#' @param function_temp_k Number of radial basis functions used to construct the
#' nonlinear temperature effect.
#'
#' @param function_temp_amplitude Standard deviation of the nonlinear
#' temperature function after rescaling.
#'
#' @param dust_alpha_rho Correlation between installation baseline efficiency
#' and installation-specific dust effects.
#'
#' @param sd_installation_dust Standard deviation of installation-specific dust
#' effects.
#'
#' @param sd_installation_alpha Standard deviation of installation baseline
#' efficiencies.
#'
#' @param baseline_installation_mean Population-average baseline log-odds of
#' efficiency.
#'
#' @param baseline_installation_dust Population-average dust effect on
#' efficiency.
#'
#' @param kappa_efficiency Precision parameter of the Beta distribution.
#' Larger values generate less variable efficiency measurements.
#'
#' @param temp_weights Numeric vector containing the spline coefficients used to
#' construct the nonlinear temperature function. Supplying different weights
#' changes the true shape of the temperature-efficiency relationship.
#'
#' @return
#' A data frame with one row per installation-day measurement containing:
#'
#' * `instalations_id` — installation identifier.
#' * `temperature` — daily air temperature (°C).
#' * `dust_index` — standardized dust accumulation index.
#' * `efficiency` — simulated solar panel efficiency (0–1).
#'
#' @examples
#' sim_data <- simulate_solar_efficiency_vs_temperature()
#'
#' plot(
#'   sim_data$temperature,
#'   sim_data$efficiency,
#'   pch = 16
#' )
#'
#' @export
simulate_solar_efficiency_vs_temperature <- function(
  num_installations = 30,
  num_days          = 10,
  mean_temperature  = 25,
  sd_temperature    = 10,
  function_temp_k   = 15,
  function_temp_amplitude = 1,
  dust_alpha_rho          = -0.5,
  sd_installation_dust    = 0.3,
  sd_installation_alpha   = 0.2,
  baseline_installation_mean = 0.5,
  baseline_installation_dust = -0.3,
  kappa_efficiency           = 30,
  temp_weights = c(0,1,3,6,9,11,12,11,9,6,3,1,0,0,0)
){

  # Simulate the solar intallations and the days 
  idx <- make_nested_ids(levels = list(instalations = num_installations, num_days = 10))

  # Simulate Diffrent Temperatures and Dust Index
  temperature <- rnorm(n = nrow(idx), mean = mean_temperature, sd = sd_temperature)
  dust_index  <- runif(n = nrow(idx), min = 0, max = 1)

  # Make function temp using Radial gaussian bumps
  f_t <- make_radial_basis(x = temperature, knots = function_temp_k, length_scale = 10)
  f_t <- simulate_spline_effect(
    B         = f_t, 
    amplitude = function_temp_amplitude,
    weights   = temp_weights
  )
  
  # Make Correlated hierarcle intercept and dust index coef
  R_cor      <- matrix(c(1, dust_alpha_rho, dust_alpha_rho, 1), nrow = 2, ncol = 2)
  cor_effect <- make_correlated_effects(
    n                  = num_installations, 
    correlation_matrix = R_cor, 
    means              = c(0, 0), 
    sds = c(sd_installation_dust, sd_installation_alpha) 
  )
  beta_dust_index <- baseline_installation_dust + cor_effect[, 1] 
  alpha_j         <- baseline_installation_mean + cor_effect[, 2]

  # Make a Linear Predictor 
  mu_i <- plogis(
    alpha_j[idx$instalations_id] +
      beta_dust_index[idx$instalations_id] * dust_index + f_t$fitted
  )
  # Sample the efficiency from Beta distribution 
  efficiency <- rbeta(
    n      = nrow(idx),
    shape1 = mu_i * kappa_efficiency,
    shape2 = (1 - mu_i) * kappa_efficiency
  )

  # Combine in one Dataframe
  sim_data <- data.frame(
    instalations_id = idx$instalations_id,
    dust_index      = dust_index,
    temperature     = temperature,
    efficiency      = efficiency
  )
}