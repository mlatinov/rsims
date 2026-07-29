#' Simulate Ski Patrol Safety Briefing Data
#'
#' @description
#' ## Scenario
#'
#' Each ski resort runs avalanche-safety briefings. The goal is to estimate
#' the effect of briefing attendance (measured as hours attended) on whether a
#' skier experiences an incident.
#'
#' Briefings are randomly assigned, meaning there is no confounding between
#' briefing attendance and incident risk. However, resorts are not identical:
#' some resorts have higher baseline incident rates, and some resorts have
#' more effective safety briefings than others.
#'
#' The challenge is to recover the average effect of briefing attendance while
#' accounting for variation between resorts. In particular, resorts with higher
#' baseline incident risk are assumed to benefit more from longer briefings.
#'
#' This creates a hierarchical causal problem where the treatment effect varies
#' across groups:
#'
#' \deqn{
#' \text{Resort Risk}_j \rightarrow
#' \text{Briefing Effect}_j
#' }
#'
#' with a negative correlation between baseline risk and briefing effectiveness.
#' Resorts with larger baseline risks tend to have stronger protective effects
#' from briefings.
#'
#' ## Data Generating Model
#'
#' For resort \eqn{j} and skier \eqn{i}:
#'
#' \deqn{
#' Y_{ij} \sim Bernoulli(p_{ij})
#' }
#'
#' where:
#'
#' \deqn{
#' logit(p_{ij}) =
#' \alpha_j +
#' \beta_j BriefingHours_{ij}
#' }
#'
#' Resort-specific intercepts and slopes are sampled jointly:
#'
#' \deqn{
#' \begin{pmatrix}
#' \alpha_j\\
#' \beta_j
#' \end{pmatrix}
#' \sim
#' MVN
#' \left(
#' \begin{pmatrix}
#' \alpha_0\\
#' \beta_0
#' \end{pmatrix},
#' \Sigma
#' \right)
#' }
#'
#' The covariance matrix is constructed as:
#'
#' \deqn{
#' \Sigma = D R D
#' }
#'
#' where \eqn{R} controls the correlation between baseline incident risk and briefing effectiveness.
#' 
#' @param num_resorts Number of ski resorts to simulate.
#' @param briefings_per_resorts Number of safety briefings performed per resort.
#' @param briefings_time_constrains Vector of minimum and maximum briefing lenghts in hours
#' @param mean_resorts_prob Average baseline log-odds of an incident across resorts.
#' @param mean_briefing_effect Average effect of briefing duration on incident risk.
#' @param sigma_briefing_effect Standard deviation of the resort-specific random effects.
#' @param rho Correlation between resort baseline risk and briefing effectiveness.
#'
#' @return A data.frame containing:
#' \itemize{
#'   \item resort_id: Identifier of the ski resort.
#'   \item briefings_hours: Duration of the safety briefing in hours.
#'   \item incident: Binary indicator of whether an incident occurred.
#' }
#'
#' @export
simulate_ski_patrol <- function(
  num_resorts = 30,
  briefings_time_constrains = c(0.30, 2),
  briefings_per_resorts = 2,
  mean_resorts_prob = 0.003,
  mean_briefing_effect = -2,
  sigma_briefing_effect = 0.5,
  rho = -0.6
){
  ## Simulate N Resorts that run safety briefings
  n <- num_resorts * briefings_per_resorts
  resort_id <- rep(seq_len(num_resorts), each = briefings_per_resorts)

  # Simulate the Briefing Length in Hours
  briefings_h <- runif(n, min = min(briefings_time_constrains), max = max(briefings_time_constrains))

  ## Correlated Varying Effects between the Baseline and the Beta Coefficient for Briefings
  R_corr <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  sigmas <- diag(c(sigma_briefing_effect, sigma_briefing_effect))
  Sigma <- sigmas %*% R_corr %*% sigmas

  # Sample from Multivariate Normal and Correlate the Coefficients
  u <- MASS::mvrnorm(n = num_resorts, mu = c(0, 0), Sigma = Sigma)
  alpha_j <- mean_resorts_prob + u[, 1]
  beta_j  <- mean_briefing_effect + u[, 2]

  # Make the Linear Predictor and Sample from Binomial Distribution
  eta      <- plogis(alpha_j[resort_id] + beta_j[resort_id] * briefings_h)
  incident <- rbinom(n, size = 1, prob = eta)

  # Combine into one Dataset
  sim_data <- data.frame(
    resort_id = resort_id,
    briefings_hours = briefings_h,
    incident = incident
  )
}