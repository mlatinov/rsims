#' Simulate Job Satisfaction Survey Data
#'
#' @description
#' Simulates a simple one-factor Confirmatory Factor Analysis (CFA) model for
#' employee job satisfaction.
#'
#' A latent job satisfaction variable is generated for every employee and is
#' measured indirectly through multiple continuous survey items. Each indicator
#' has its own intercept, loading, and measurement error.
#'
#' The data-generating process is
#'
#' \deqn{
#' \eta_i \sim Normal(0,1),
#' }
#'
#' and each observed survey item follows
#'
#' \deqn{
#' y_{ij}
#' =
#' \nu_j
#' +
#' \lambda_j \eta_i
#' +
#' \varepsilon_{ij},
#' }
#'
#' where
#'
#' \deqn{
#' \varepsilon_{ij}
#' \sim
#' Normal(0,\sigma_j).
#' }
#'
#' The latent variance is fixed to one, corresponding to the standardised
#' identification strategy commonly used in Confirmatory Factor Analysis.
#'
#' This simulation is useful for teaching and benchmarking CFA models because
#' the latent variable is known exactly before measurement error is added.
#'
#' @param employees Integer.
#' Number of employees to simulate.
#'
#' @param v Numeric vector.
#' Indicator intercepts (\eqn{\nu}).
#'
#' @param measurement_noise Numeric vector.
#' Residual standard deviations for each indicator.
#'
#' @param lambda Numeric vector.
#' Factor loadings relating the latent job satisfaction variable to each
#' indicator.
#'
#' @param indicator_names Character vector.
#' Names assigned to the simulated survey items.
#'
#' @return
#' A data frame containing one row per employee with the simulated questionnaire
#' responses.
#'
#' @details
#' The simulation corresponds to a single-factor Confirmatory Factor Analysis
#' model:
#'
#' \itemize{
#' \item One latent variable (job satisfaction).
#' \item Continuous indicators.
#' \item Conditional independence of indicators given the latent variable.
#' \item Standardised latent scaling (\eqn{Var(\eta)=1}).
#' }
#'
#' Because the latent variable is generated before the indicators, the true
#' measurement model is known exactly and can be used to validate CFA software
#' such as **lavaan**, **blavaan**, or **Stan** implementations.
#'
#' @examples
#' sim <- simulate_latent_job_satisfaction()
#'
#' head(sim)
#'
#' cor(sim)
#'
#' @references
#' Brown, T. A. (2015).
#' *Confirmatory Factor Analysis for Applied Research* (2nd ed.).
#'
#' Bollen, K. A. (1989).
#' *Structural Equations with Latent Variables*.
#'
#' @export
simulate_latent_job_satisfaction <- function(
  employees = 1200,
  v                 = c(4.0, 3.5, 4.2),
  measurement_noise = c(0.5, 0.6, 0.45),
  lambda            = c(0.9, 0.7, 1.1),
  indicator_names   = c("Q1", "Q2", "Q3")
){

  # Simulate job satisfactions scores from the indicators 
  latent_job_satisfaction <- simulate_latent(n = employees, mean = 0, sd = 1)
  
  # Simulate the measurments back  
  measurments_job_satisfaction <- simulate_indicators(
    eta      = latent_job_satisfaction$data$eta,
    loadings = lambda,
    sigma    = measurement_noise,
    intercepts = v,
    names      = indicator_names
  )

  # Return the recorved data from the measurments 
  sim_data <- data.frame(measurments_job_satisfaction$data)
}