#' Simulate Body Temperature Dynamics Using a First-Order ODE
#'
#' Simulates repeated body temperature measurements for multiple subjects using
#' a mechanistic ordinary differential equation (ODE). Each subject has a
#' subject-specific heat production rate and heat loss rate, producing
#' heterogeneous temperature trajectories over time.
#'
#' The underlying process follows
#'
#' \deqn{
#' \frac{dT}{dt} = a_j - b_jT
#' }
#'
#' where \eqn{T} denotes body temperature, \eqn{a_j} is the subject-specific
#' heat inflow rate, and \eqn{b_j} is the subject-specific heat loss rate.
#'
#' The analytical solution is
#'
#' \deqn{
#' T_j(t)
#' =
#' \frac{a_j}{b_j}
#' +
#' \left(
#' T_{0,j}
#' -
#' \frac{a_j}{b_j}
#' \right)
#' e^{-b_jt}.
#' }
#'
#' Heat inflow depends on ambient temperature:
#'
#' \deqn{
#' a_j
#' =
#' \exp\left(
#' \mu_a
#' +
#' \beta_{ambient}A_j
#' +
#' \epsilon_j
#' \right)
#' }
#'
#' while heat loss depends on physical activity:
#'
#' \deqn{
#' b_j
#' =
#' \exp\left(
#' \mu_b
#' +
#' \beta_{activity}L_j
#' +
#' \eta_j
#' \right)
#' }
#'
#' Gaussian measurement error is added to the analytical solution to produce
#' observed body temperatures.
#'
#' @param num_subjects Number of subjects.
#' @param num_hours Number of hourly observations per subject.
#' @param mean_ambient_t Mean ambient temperature.
#' @param sd_ambient_t Standard deviation of ambient temperature.
#' @param mean_activity_l Mean activity level.
#' @param sd_activity_l Standard deviation of activity level.
#' @param mean_initial_point Mean initial body temperature.
#' @param sd_initial_point Standard deviation of initial body temperature.
#' @param mean_heat_inflow_rate Mean log heat inflow rate.
#' @param sd_subjects_heat_inflow_rate Between-subject standard deviation of heat inflow.
#' @param beta_ambient_t Effect of ambient temperature on heat inflow.
#' @param mean_outflow_rate Mean log heat loss rate.
#' @param sd_subjects_outflow_rate Between-subject standard deviation of heat loss.
#' @param beta_activity Effect of activity level on heat loss.
#' @param body_temp_sd Observation error standard deviation.
#'
#' @return A data frame with one row per subject-hour containing:
#' \describe{
#'   \item{hour}{Hour of observation.}
#'   \item{subject}{Subject identifier.}
#'   \item{ambient_temp}{Ambient temperature for the subject.}
#'   \item{activity_level}{Subject activity level.}
#'   \item{body_temp}{Simulated observed body temperature.}
#' }
#'
#' @details
#' This simulator illustrates a simple mechanistic model in which body
#' temperature approaches a subject-specific equilibrium determined by the
#' balance between heat production and heat loss. The model demonstrates how
#' ordinary differential equations can be used to generate longitudinal data
#' with interpretable physiological parameters.
#'
#' @examples
#' data <- simulate_ode_body_temperature()
#'
#' head(data)
#'
#' @export
simulate_ode_body_temperature <- function(
    num_subjects       = 40,
    num_hours          = 48,
    mean_ambient_t = 20,               sd_ambient_t   = 4,
    mean_activity_l  = 5,              sd_activity_l    = 1.5,
    mean_initial_point = 36.8,         sd_initial_point = 0.5,
    mean_outflow_rate = log(0.15),   
    sd_subjects_outflow_rate = 0.08,
    beta_activity = -0.03,
    mean_heat_inflow_rate = log(0.15) + log(37), 
    sd_subjects_heat_inflow_rate = 0.08,
    beta_ambient_t = 0.015,
    body_temp_sd = 0.15
){
  # Simulate Subjects measurment of body temperature every hour to max hours
  n   <- num_subjects * num_hours 
  ids <- make_nested_ids(levels = list(subject = num_subjects, hour = num_hours)) 

  # Covariates 
  ambient_temp   <- rnorm(num_subjects, mean = mean_ambient_t, sd = sd_ambient_t)
  activity_level <- rnorm(num_subjects, mean = mean_activity_l, sd = sd_activity_l)

  # Center the covariates 
  ambient_c  <- ambient_temp   - mean_ambient_t
  activity_c <- activity_level - mean_activity_l

  # Paramters 
  a  <- exp(mean_heat_inflow_rate + beta_ambient_t * ambient_c + sd_subjects_heat_inflow_rate * rnorm(num_subjects, 0, 1))
  b  <- exp(mean_outflow_rate     + beta_activity  * activity_c + sd_subjects_outflow_rate * rnorm(num_subjects, 0 , 1))
  xo <- rnorm(num_subjects, mean = mean_initial_point, sd = sd_initial_point)

  # ODE Predictor 
  mu <- (
    (a[ids$subject_id] / b[ids$subject_id]) + (xo[ids$subject_id] - (a[ids$subject_id] / b[ids$subject_id]))
    * exp(-b[ids$subject_id] * ids$hour_id)
  )
  # Sample from normal distribution 
  body_temp <- rnorm(n, mean = mu, sd = body_temp_sd)

  # Return a dataframe 
  data.frame(
    hour = ids$hour_id,
    subject = ids$subject_id,
    ambient_temp = ambient_temp[ids$subject_id],
    activity_level = activity_level[ids$subject_id],
    body_temp = body_temp
  )
}