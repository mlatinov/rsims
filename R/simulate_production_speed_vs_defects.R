#' Simulate Machine Speed Effects on Manufacturing Defects
#'
#' @description
#' Simulates a hierarchical manufacturing process where multiple factories
#' operate machines at different speeds across repeated shifts. The number of
#' defects is generated from a Poisson distribution with a nonlinear effect of
#' operating speed.
#'
#' Each factory has its own baseline defect rate through a hierarchical
#' intercept:
#'
#' \deqn{
#' \alpha_j = \mu_\alpha + \sigma_\alpha z_j,
#' \qquad z_j \sim Normal(0,1)
#' }
#'
#' The nonlinear relationship between machine speed and defects is generated
#' using a radial basis spline representation:
#'
#' \deqn{
#' f(speed_i)=B(speed_i)w
#' }
#'
#' where the basis matrix is created from radial Gaussian basis functions and
#' the spline weights control the shape of the speed-defect relationship.
#'
#' The expected number of defects follows a Poisson log-link model:
#'
#' \deqn{
#' \lambda_i =
#' \exp(\alpha_j + f(speed_i))
#' }
#'
#' and observed defects are sampled as:
#'
#' \deqn{
#' defects_i \sim Poisson(\lambda_i)
#' }
#'
#' This simulation is useful for testing hierarchical count models,
#' nonlinear regression models, Poisson regression, and spline-based models.
#'
#' @param num_factories Number of independent factories to simulate.
#'
#' @param num_shifts_per_factory Number of production shifts observed for each
#' factory.
#'
#' @param mean_factories_defects Mean log-scale baseline defect rate across
#' factories.
#'
#' @param sd_factories_defects Standard deviation of factory-level baseline
#' defect variation.
#'
#' @param speed_weights Numeric vector defining the spline weights controlling
#' the nonlinear effect of machine speed. Larger values create stronger
#' increases in defects at higher speeds.
#'
#' @param speed_sigma_group Standard deviation controlling factory-level
#' variation in the nonlinear speed effect.
#'
#' @return
#' A data frame containing:
#'
#' \itemize{
#'   \item \code{factory_id}: Identifier for the factory.
#'   \item \code{speed}: Machine operating speed during the shift.
#'   \item \code{defects_count}: Number of defects observed during the shift.
#' }
#'
#' @examples
#' set.seed(123)
#'
#' sim_data <- simulate_production_speed_vs_defects()
#'
#' plot(
#'   sim_data$speed,
#'   sim_data$defects_count,
#'   xlab = "Machine Speed",
#'   ylab = "Number of Defects"
#' )
#'
#' @export
simulate_production_speed_vs_defects <- function(
  num_factories = 5,
  num_shifts_per_factory = 60,
  mean_factories_defects = log(5),
  sd_factories_defects   = 0.5,
  speed_weights = c(0,0,0.2,0.8,2,5,10,18),
  speed_sigma_group = 0.5
){
  # Simulate the J factories and N shifts per factory 
  idx <- make_nested_ids(levels = c(factory = num_factories, shift = num_shifts_per_factory))

  # Simulate operating machine speed 
  speed <- runif(n = nrow(idx), min = 10, max = 60)

  # Defects paramters 
  alpha_j <- mean_factories_defects + sd_factories_defects * rnorm(n = num_factories, mean = 0, sd = 1)
  B       <- make_radial_basis(x = speed, knots = 8, length_scale = 12) 
  f_speed <- simulate_hierarchical_spline_effect(
    amplitude = 1,
    B = B,
    group = idx$factory_id,
    weights     = speed_weights,
    sigma_group = speed_sigma_group
  )
  # Create a linear predictor on a exp scale and sample from Possion Distribution 
  lambda_i  <- exp(alpha_j[idx$factory_id] + f_speed$fitted)
  defects_count <- rpois(n = nrow(idx), lambda = lambda_i) 

  # Combine in one dataframe 
  sim_data <- data.frame(
    factory_id = idx$factory_id,
    speed      = speed,
    defects_count = defects_count
  )
}