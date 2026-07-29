#' Simulate Athlete Fitness Data
#'
#' @description
#' Simulates a Multiple Indicators Multiple Causes (MIMIC) model for athlete
#' fitness.
#'
#' Athlete fitness is represented as a latent variable predicted by observed
#' covariates (training hours, age, and sleep duration). The latent fitness
#' variable is then measured indirectly through several continuous fitness
#' assessments.
#'
#' The structural model is
#'
#' \deqn{
#' \eta_i =
#' \beta_{age}Age_i +
#' \beta_{train}Training_i +
#' \beta_{sleep}Sleep_i +
#' \zeta_i,
#' }
#'
#' where
#'
#' \deqn{
#' \zeta_i \sim Normal(0,\sigma_{\zeta}).
#' }
#'
#' Each observed indicator follows
#'
#' \deqn{
#' y_{ij}
#' =
#' \nu_j
#' +
#' \lambda_j\eta_i
#' +
#' \varepsilon_{ij},
#' }
#'
#' with
#'
#' \deqn{
#' \varepsilon_{ij}
#' \sim
#' Normal(0,\sigma_j).
#' }
#'
#' This corresponds to a MIMIC model in which the observed covariates influence
#' the fitness assessments only through the latent fitness variable.
#'
#' @param n Integer.
#' Number of athletes to simulate.
#'
#' @param mean_train_h Numeric.
#' Mean weekly training hours.
#'
#' @param age_range Numeric vector of length two.
#' Minimum and maximum athlete age.
#'
#' @param mean_sleed_h Numeric.
#' Mean nightly sleep duration in hours.
#'
#' @param beta_age Numeric.
#' Effect of age on latent fitness.
#'
#' @param beta_train Numeric.
#' Effect of training hours on latent fitness.
#'
#' @param beta_sleep Numeric.
#' Effect of sleep duration on latent fitness.
#'
#' @param lambda Numeric vector.
#' Factor loadings for the fitness indicators.
#'
#' @param measurment_noise Numeric vector.
#' Residual standard deviations of the indicators.
#'
#' @param intercepts Numeric vector.
#' Indicator intercepts.
#'
#' @param indicator_names Character vector.
#' Names assigned to the simulated fitness tests.
#'
#' @return
#' A data frame containing
#'
#' \itemize{
#' \item Three continuous fitness indicators.
#' \item Athlete age.
#' \item Weekly training hours.
#' \item Average sleep duration.
#' }
#'
#' @details
#' This simulation represents a classical MIMIC model:
#'
#' \itemize{
#' \item Observed covariates predict a latent variable.
#' \item The latent variable is measured by multiple continuous indicators.
#' \item Covariates have no direct effects on the indicators.
#' \item Identification uses the marker-variable approach because the latent
#' variance is determined by the structural regression.
#' }
#'
#' By default, the indicators can be interpreted as measurements such as
#' VO2 max, shuttle-run performance, and sprint performance.
#'
#' This simulation is useful for validating Bayesian and frequentist MIMIC
#' models, Structural Equation Models (SEM), and latent regression models.
#'
#' @examples
#' sim <- simulate_athlete_fitness()
#'
#' head(sim)
#'
#' cor(sim)
#'
#' @references
#' Bollen, K. A. (1989).
#' *Structural Equations with Latent Variables*.
#'
#' Brown, T. A. (2015).
#' *Confirmatory Factor Analysis for Applied Research* (2nd ed.).
#'
#' Kline, R. B. (2023).
#' *Principles and Practice of Structural Equation Modeling* (5th ed.).
#'
#' @export
simulate_athlete_fitness <- function(
  n = 900,
  # Covariates 
  mean_train_h = 3,
  age_range    = c(20, 40),
  mean_sleep_h = 6,
  beta_age = -0.03,
  beta_train = 0.05,
  beta_sleep = 0.12, 
  # Indicators 
  lambda            = c(6.0, 2.0, 0.35),
  measurement_noise = c(3.0, 1.2, 0.25),
  intercepts        = c(45, 12, 2.5),
  indicator_names = c("vo2_max","shuttle_run","sprint_index")
){
  # Simulate the covariates 
  age     <- runif(n, min = age_range[1], max = age_range[2])
  train_h <- rnorm(n, mean = mean_train_h, sd = 1)
  sleep_h <- rnorm(n, mean = mean_sleep_h, sd = 1)

  # Latent Regreesion towards Eta fitness performance 
  fitness_performance <- simulate_latent_regression(
    X          = data.frame(age, train_h, sleep_h),
    betas      =  c(beta_age, beta_train, beta_sleep),
    sigma_zeta = 1
  )

  # Simulate indicators 
  indicators <- simulate_indicators(
    eta      = fitness_performance$data$eta,
    loadings = lambda,
    intercepts = intercepts,
    sigma = measurement_noise,
    names = indicator_names
  )

  # Combine the indicators with the covariates and return it 
  sim_data <- cbind(indicators$data$Y, data.frame(age, train_h, sleep_h))
}