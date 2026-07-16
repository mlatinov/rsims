#' Simulate Fuel Efficiency as a Nonlinear Function of Speed
#'
#' @description
#' Simulates vehicle fuel efficiency measurements collected across different
#' driving speeds. The relationship between speed and fuel efficiency is
#' generated using a smooth P-spline nonlinear function rather than a fixed
#' parametric equation.
#'
#' The simulation represents a scenario where the true relationship between
#' speed and fuel efficiency is unknown and must be learned from the data.
#'
#' The generated model is:
#'
#' \deqn{
#' FuelEfficiency_i \sim Normal(\mu_i,\sigma)
#' }
#'
#' where the expected fuel efficiency is:
#'
#' \deqn{
#' \mu_i = \alpha + f(speed_i)
#' }
#'
#' and the nonlinear speed effect is generated as:
#'
#' \deqn{
#' f(speed_i)=\sum_{k=1}^{K}B_k(speed_i)w_k
#' }
#'
#' where:
#'
#' \itemize{
#'   \item \eqn{B_k(speed_i)} are precomputed spline basis functions,
#'   \item \eqn{w_k} are spline weights generated through a random walk,
#'   \item \eqn{\alpha} is the baseline fuel efficiency.
#' }
#'
#' The spline weights follow a first-order P-spline random walk:
#'
#' \deqn{
#' w_1=\tau z_1
#' }
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
#' The parameter \eqn{\tau} controls the smoothness of the underlying
#' relationship between speed and fuel efficiency.
#'
#' @param num_test_drives Number of simulated vehicle tests.
#'
#' @param speed_ranges Numeric vector of length two defining the minimum and
#' maximum possible driving speeds.
#'
#' @param tau Random walk step size controlling spline smoothness.
#' Smaller values generate smoother speed-efficiency relationships.
#'
#' @param baseline_fuel_efficiency Baseline fuel efficiency value when the
#' nonlinear speed contribution is centered around zero.
#'
#' @param sd_fuel_efficiency Residual standard deviation of fuel efficiency
#' measurements.
#'
#' @param B_df Number of spline basis functions used to construct the nonlinear
#' speed effect.
#'
#' @return
#' A data frame containing:
#'
#' \describe{
#'   \item{speed}{The simulated driving speed.}
#'   \item{fuel_efficiency}{Observed fuel efficiency measurements generated
#'   from the nonlinear model.}
#' }
#'
#' @examples
#'
#' data <- simulate_fuel_efficiency_vs_speed(
#'   num_test_drives = 500
#' )
#'
#' plot(
#'   data$speed,
#'   data$fuel_efficiency
#' )
#'
#' @export
simulate_fuel_efficiency_vs_speed <- function(
  num_test_drives = 500,
  speed_ranges    = c(20, 160),
  tau             = 0.15,
  baseline_fuel_efficiency = 20,
  sd_fuel_efficiency       = 1,
  B_df = 8
){
  # Simulate Diffrent speed test across the number of test drives 
  speed <- runif(num_test_drives, min = speed_ranges[1], max = speed_ranges[2])

  # Create a Basis P Spline function of speed 
  B  <- splines::bs(x = speed, df = B_df)
  k  <- ncol(B)
  zk <- rnorm(n = k, mean = 0, sd = 1) 
  w  <- numeric(k)
  w[1] <- tau * zk[1]
  for(i in 2:k){
    w[i] <- w[i-1] + tau * zk[i]
  }
  # Create the f out of speed 
  fs <- B %*% w
  # rescale the nonlinear effect
  fs <- fs / sd(fs) * 5

  ## Sample Fuel Efficiency from Normal Distribution 
  # Linear predictor 
  mu_i             <- baseline_fuel_efficiency + fs 
  fuel_efficiency  <- rnorm(n = num_test_drives, mean = mu_i, sd = sd_fuel_efficiency)

  # Combine in one simulated dataset 
  sim_data <- data.frame(
    speed = speed,
    fuel_efficiency = fuel_efficiency
  )
}