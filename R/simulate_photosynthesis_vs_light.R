#' Simulate Photosynthesis Response to Light Intensity
#'
#' @description
#' Simulates leaf-level photosynthesis measurements as a nonlinear function of
#' incident light intensity (PAR). The relationship between light availability
#' and CO2 assimilation is generated using a smooth P-spline function.
#'
#' This represents a biological scenario where photosynthetic response does not
#' follow a simple linear relationship. At different light levels, leaves may
#' show changing rates of carbon assimilation due to saturation and nonlinear
#' physiological responses.
#'
#' The data generating model is:
#'
#' \deqn{
#' CO2_i \sim Normal(\mu_i,\sigma)
#' }
#'
#' with:
#'
#' \deqn{
#' \mu_i = \alpha + f(light_i)
#' }
#'
#' where the nonlinear light response is:
#'
#' \deqn{
#' f(light_i)=\sum_{k=1}^{K}B_k(light_i)w_k
#' }
#'
#' The spline coefficients follow a P-spline random walk prior:
#'
#' \deqn{
#' w_k=w_{k-1}+\tau z_k
#' }
#'
#' where:
#'
#' \deqn{
#' z_k \sim Normal(0,1)
#' }
#'
#' The raw spline effect is standardized and scaled by
#' \code{spline_amplitude} to control the biological magnitude of the nonlinear
#' response.
#'
#' @param num_leaf_measurements Number of simulated leaf measurements.
#'
#' @param light_par_range Numeric vector of length two defining the minimum and
#' maximum PAR light intensity values.
#'
#' @param tau_light Controls the smoothness of the P-spline random walk.
#' Smaller values generate smoother light response curves.
#'
#' @param spline_amplitude Controls the magnitude of the nonlinear light effect.
#'
#' @param baseline_co2_assimilation_rate Baseline CO2 assimilation rate around
#' which the nonlinear light effect varies.
#'
#' @param sd_co2_assimilation_rate Measurement noise standard deviation.
#'
#' @return
#' A data frame containing:
#'
#' \describe{
#'   \item{light_par}{Simulated PAR light intensity.}
#'   \item{co2_assimilation_rate}{Observed CO2 assimilation measurements.}
#'   \item{true_mu}{True expected assimilation rate before observation noise.}
#'   \item{true_light_effect}{True nonlinear contribution from light.}
#' }
#'
#' @examples
#'
#' data <- simulate_photosynthesis_vs_light()
#'
#' plot(
#'   data$light_par,
#'   data$co2_assimilation_rate
#' )
#'
#' @export
simulate_photosynthesis_vs_light <- function(
  num_leaf_measurements = 500,
  light_par_range = c(100, 2000),
  tau_light = 0.15,
  spline_amplitude = 8,
  baseline_co2_assimilation_rate = 10,
  sd_co2_assimilation_rate = 1
){
  # Simulate Light Par values for every leaf
  light <- runif(n = num_leaf_measurments, min = light_par_range[1], max = light_par_range[2])

  # Build a Non linear P-Spline function of light 
  f_light <- simulate_p_spline(x = light, df = spline_amplitude, tau = tau_light)

  # Sample the CO2 assimilation rate from Normal Distribution 
  mu_i                  <- baseline_co2_assimilation_rate + f_light$fitted
  co2_assimilation_rate <- pmax(rnorm(n = num_leaf_measurments, mean = mu_i, sd = sd_co2_assimilation_rate),0)

  # Combine in one dataset 
  sim_data <- data.frame(
    light_par = light,
    co2_assimilation_rate = co2_assimilation_rate
  )
}

