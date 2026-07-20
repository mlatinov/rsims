#' Simulate Plant Vigour and Stress Tolerance Measurements
#'
#' @description
#' Simulates a two-factor Confirmatory Factor Analysis (CFA) model describing
#' plant vigour and stress tolerance.
#'
#' Two correlated latent variables are generated for each plant and are measured
#' through separate groups of continuous indicators. The latent variables are
#' correlated but neither causes the other.
#'
#' The latent variables follow
#'
#' \deqn{
#' (\eta_{vig},\eta_{tol})
#' \sim
#' MVN(\mathbf{0}, \mathbf{R}),
#' }
#'
#' where both latent variances are fixed to one and
#' \eqn{\mathbf{R}} contains the latent correlation.
#'
#' Each observed indicator follows
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
#' with
#'
#' \deqn{
#' \varepsilon_{ij}
#' \sim
#' Normal(0,\sigma_j).
#' }
#'
#' Indicators are conditionally independent given their corresponding latent
#' variable.
#'
#' @param n Integer.
#' Number of plants to simulate.
#'
#' @param v Numeric vector.
#' Indicator intercepts for the plant vigour latent variable.
#'
#' @param lambda_v Numeric vector.
#' Factor loadings for the vigour indicators.
#'
#' @param measurement_noise_v Numeric vector.
#' Residual standard deviations of the vigour indicators.
#'
#' @param names_v Character vector.
#' Names assigned to the vigour indicators.
#'
#' @param t Numeric vector.
#' Indicator intercepts for the stress tolerance latent variable.
#'
#' @param lambda_t Numeric vector.
#' Factor loadings for the stress tolerance indicators.
#'
#' @param measurement_noise_t Numeric vector.
#' Residual standard deviations of the stress tolerance indicators.
#'
#' @param names_t Character vector.
#' Names assigned to the stress tolerance indicators.
#'
#' @param rho Numeric.
#' Correlation between the latent plant vigour and stress tolerance variables.
#'
#' @return
#' A data frame containing six simulated continuous measurements:
#'
#' \itemize{
#' \item Three indicators measuring plant vigour.
#' \item Three indicators measuring stress tolerance.
#' }
#'
#' @details
#' The simulated model represents a classical two-factor Confirmatory Factor
#' Analysis (CFA):
#'
#' \itemize{
#' \item Two correlated latent variables.
#' \item Three continuous indicators per latent factor.
#' \item Conditional independence of indicators given their latent variable.
#' \item Standardised latent variances (\eqn{Var(\eta)=1}).
#' }
#'
#' The default indicators correspond to common plant physiology measurements.
#' The vigour factor is measured by shoot height, leaf count and stem diameter,
#' while the stress tolerance factor is measured by Fv/Fm ratio, proline
#' concentration and relative water content.
#'
#' This simulation is useful for validating multi-factor CFA, Bayesian latent
#' variable models, Structural Equation Models (SEM), and multivariate
#' measurement models.
#'
#' @examples
#' sim <- simulate_plant_vigour_stress_tolerance()
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
simulate_plant_vigour_stress_tolerance <- function(
  n = 2000,
  # Vigour indicators
  v        = c(25, 12, 3.0),
  lambda_v = c(4.0, 2.5, 0.6),
  measurement_noise_v = c(2.0, 1.5, 0.3),
  names_v = c("V1", "V2", "V3"),
  # Tolerance indicators
  t        = c(0.75, 18, 55),
  lambda_t = c(0.08, 3.0, 7.0),
  measurement_noise_t = c(0.05, 2.0, 5.0),
  names_t = c("T1", "T2", "T3"),
  # Corrlation Latent 
  rho = 0.6
){

  # Simulate Vigour and Tolerance latent variables with corrleated
  cor_matrix <- matrix(c(1 ,rho, rho, 1), ncol = 2) 
  latents    <- simulate_correlated_latents(
    n = n, 
    sds = c(1, 1),
    Rho = cor_matrix,
    names = c("vigour","stress_tolerance")
  )

  # Convert to loadins to a matrix format 
  lambda <- matrix(0, nrow = 2, ncol = 6)
  lambda[1, 1:3] <- lambda_v
  lambda[2, 4:6] <- lambda_t

  # Simulate back the measurments indicators 
  indicators <- simulate_indicators(
    eta      = latents$data$Eta,
    loadings = lambda,
    intercepts = c(v, t),
    sigma      = c(measurement_noise_v, measurement_noise_t),
    names      = c(names_v, names_t)
  )
  # Convert to a dataframe and return 
  as.data.frame(x = indicators$data$Y)
}