#' Generate Independent Random Effects
#'
#' @description
#' Generates normally distributed group-level effects.
#'
#' The generated values follow:
#'
#' \deqn{
#' \theta_j \sim Normal(\mu,\sigma)
#' }
#'
#' where each group receives its own random deviation.
#'
#' These effects are commonly used for varying intercepts or varying slopes
#' in hierarchical models.
#'
#' @param n Number of groups.
#' @param mean Population mean of the random effect.
#' @param sd Standard deviation of the random effect.
#'
#' @return
#' A numeric vector containing one random effect per group.
#'
#' @examples
#'
#' hospital_effects <- make_random_effects(
#'   n = 20,
#'   mean = 0,
#'   sd = 2
#' )
#'
#' @keywords internal
make_random_effects <- function(
  n,
  mean = 0,
  sd = 1
){

  mean + sd * rnorm(
    n = n,
    mean = 0,
    sd = 1
  )

}



#' Generate Correlated Random Effects
#'
#' @description
#' Generates multivariate normally distributed random effects with a user
#' specified covariance structure.
#'
#' The generated effects follow:
#'
#' \deqn{
#' \mathbf{u}_j \sim MVN(\boldsymbol{\mu},\Sigma)
#' }
#'
#' where:
#'
#' \deqn{
#' \Sigma = D R D
#' }
#'
#' with \eqn{D} containing the standard deviations and \eqn{R} the
#' correlation matrix.
#'
#' This is useful for simulating correlated varying intercepts and slopes.
#'
#' Examples:
#'
#' \deqn{
#' (\alpha_j,\beta_j)
#' \sim MVN
#' }
#'
#' allowing hospitals, farms, clinics, or other groups to have correlated
#' baseline values and treatment effects.
#'
#' @param n Number of groups.
#' @param means Vector of population means for each random effect.
#' @param sds Vector of standard deviations for each random effect.
#' @param correlation_matrix Correlation matrix between effects.
#'
#' @return
#' A matrix with one row per group and one column per random effect.
#'
#' @examples
#'
#' effects <- make_correlated_effects(
#'   n = 30,
#'   means = c(50, 0.5),
#'   sds = c(5, 0.2),
#'   correlation_matrix = matrix(
#'     c(1, 0.5,
#'       0.5, 1),
#'     nrow = 2
#'   )
#' )
#'
#' @keywords internal
make_correlated_effects <- function(
  n,
  means,
  sds,
  correlation_matrix
){

  if(length(means) != length(sds)){
    stop("means and sds must have the same length.")
  }

  if(nrow(correlation_matrix) != length(means)){
    stop("Correlation matrix dimensions do not match effects.")
  }


  D <- diag(sds)

  Sigma <- D %*%
    correlation_matrix %*%
    D


  MASS::mvrnorm(
    n = n,
    mu = means,
    Sigma = Sigma
  )

}