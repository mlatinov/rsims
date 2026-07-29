#' Simulate Dental Clinic Satisfaction Data
#'
#' @description
#' ## Scenario
#'
#' A network of dental clinics wants to understand which factors influence
#' patient satisfaction following an appointment.
#'
#' For every patient, you observe:
#'
#' * the number of previous visits to the clinic,
#' * whether the appointment was an emergency,
#' * the waiting time before treatment (minutes), and
#' * the patient's satisfaction score after the appointment.
#'
#' Emergency appointments are assumed to receive priority, leading to shorter
#' waiting times on average. Patients with a longer treatment history typically
#' require more complex procedures and therefore tend to wait longer.
#'
#' Clinics are not identical. Some clinics consistently receive higher
#' satisfaction scores than others, and the impact of waiting time on patient
#' satisfaction differs between clinics. Clinics with higher baseline
#' satisfaction are assumed to be less negatively affected by long waiting
#' times.
#'
#' The challenge is to estimate the effect of waiting time on patient
#' satisfaction while accounting for patient characteristics and variation
#' between clinics.
#'
#' ## Causal Structure
#'
#' \preformatted{
#'
#'                 Previous Visits ───────────────► Satisfaction
#'                       │
#'                       ▼
#'                  Waiting Time ─────────────────► Satisfaction
#'                       ▲
#'                       │
#' Emergency Status ─────┘
#'        │
#'        └──────────────────────────────────────► Satisfaction
#'
#' }
#'
#' Clinics differ in both their baseline satisfaction scores and the effect
#' of waiting time on satisfaction.
#'
#' ## Statistical Task
#'
#' Estimate the effect of waiting time on patient satisfaction while adjusting
#' for emergency appointments and previous visits. Compare a standard linear
#' regression with a hierarchical model allowing clinic-specific intercepts
#' and waiting-time effects.
#'
#' ## Data Generating Model
#'
#' Waiting times are generated as
#'
#' \deqn{
#' WaitTime_i \sim
#' Poisson(\lambda_i)
#' }
#'
#' where
#'
#' \deqn{
#' \log(\lambda_i)
#' =
#' \alpha
#' +
#' \beta_{Visit}Visits_i
#' +
#' \beta_{Emergency}Emergency_i
#' }
#'
#' Satisfaction scores are then generated as
#'
#' \deqn{
#' Satisfaction_{ij}
#' \sim
#' Normal(\mu_{ij},\sigma)
#' }
#'
#' where
#'
#' \deqn{
#' \mu_{ij}
#' =
#' \alpha_j
#' +
#' \beta_{Wait,j}WaitTime_{ij}
#' +
#' \beta_{Emergency}Emergency_{ij}
#' +
#' \beta_{Visit}Visits_{ij}
#' }
#'
#' Clinic-specific intercepts and waiting-time effects are sampled jointly:
#'
#' \deqn{
#' \begin{pmatrix}
#' \alpha_j\\
#' \beta_{Wait,j}
#' \end{pmatrix}
#' \sim
#' MVN
#' \left(
#' \begin{pmatrix}
#' \alpha_0\\
#' \beta_{Wait}
#' \end{pmatrix},
#' \Sigma
#' \right)
#' }
#'
#' where the covariance matrix is constructed as
#'
#' \deqn{
#' \Sigma = D R D
#' }
#'
#' and \eqn{R} determines the correlation between clinic baseline satisfaction
#' and the effect of waiting time.
#'
#' @param num_clinics Integer. Number of dental clinics to simulate.
#'
#' @param patients_per_clinic Integer. Number of patients generated for each
#' clinic.
#'
#' @param lambda_visits Numeric. Average number of previous visits per patient.
#'
#' @param wait_alpha Numeric. Population-average log expected waiting time
#' before accounting for previous visits and emergency status.
#'
#' @param wait_beta_visit_num Numeric. Effect of each additional previous visit
#' on the log expected waiting time. Positive values produce longer waiting
#' times for returning patients.
#'
#' @param wait_beta_emergency Numeric. Effect of emergency appointments on the
#' log expected waiting time. Negative values cause emergency patients to be
#' seen sooner.
#'
#' @param beta_emergency Numeric. Direct effect of emergency appointments on
#' patient satisfaction after accounting for waiting time.
#'
#' @param rho Numeric. Correlation between clinic-specific baseline
#' satisfaction and the clinic-specific waiting-time effect.
#'
#' @param alpha_sigma Numeric. Standard deviation of clinic-specific baseline
#' satisfaction scores.
#'
#' @param wait_time_effect_sigma Numeric. Standard deviation of the
#' clinic-specific waiting-time effects.
#'
#' @param mean_satisfaction_score Numeric. Average baseline patient
#' satisfaction across all clinics.
#'
#' @param mean_wait_time_effect Numeric. Average effect of each additional
#' minute of waiting on patient satisfaction.
#'
#' @param beta_visits Numeric. Direct effect of each previous visit on patient
#' satisfaction after adjusting for waiting time.
#'
#' @param sigma_satisfaction Numeric. Residual standard deviation of patient
#' satisfaction scores.
#'
#' @return A data frame with one row per patient containing:
#'
#' \describe{
#'   \item{clinic_id}{Unique identifier of the dental clinic.}
#'   \item{visits}{Number of previous visits by the patient.}
#'   \item{is_emergency}{Indicator for an emergency appointment (0 = no, 1 = yes).}
#'   \item{wait_time}{Waiting time before treatment (minutes).}
#'   \item{satisfaction_score}{Patient satisfaction score after the appointment.}
#' }
#'
#' @export
simulate_dental_clinics <- function(
  num_clinics = 30,
  patients_per_clinic = 40,
  lambda_visits = 2.5,
  wait_alpha = log(15),
  wait_beta_visit_num = 0.08,
  wait_beta_emergency = -0.80,
  beta_emergency = -4,
  rho = 0.5,
  alpha_sigma = 4,
  wait_time_effect_sigma = 0.20,
  mean_satisfaction_score = 80,
  mean_wait_time_effect = -0.60,
  beta_visits = -0.80,
  sigma_satisfaction = 6
){

  # Simulate J clinics each with N patients 
  n         <- num_clinics * patients_per_clinic
  clinic_id <- rep(seq_len(num_clinics), each = patients_per_clinic) 

  # Simulate Visit number 
  visit_num <- rpois(n, lambda = lambda_visits)

  # Simulate Emergency 
  emergency <- rbinom(n, size = 1, prob = 0.1)

  # Simulate the Wait time as dependent of the visit number and emergency status 
  wait_time_lambda <- exp(wait_alpha + wait_beta_visit_num * visit_num + wait_beta_emergency * emergency)
  wait_time <- rpois(n, lambda = wait_time_lambda)

  # Corrleate the effect of the baseline waiting time alphaj and the wait time 
  R_cor  <- matrix(data = c(1, rho, rho, 1),nrow = 2, ncol = 2)
  sigmas <- diag(c(alpha_sigma, wait_time_effect_sigma))
  Sigma  <- sigmas %*% R_cor %*% sigmas 

  # Sample from Multivariate Normal Distribution the alpha j and beta_wait_time j 
  u <- MASS::mvrnorm(n = num_clinics, mu = c(0, 0), Sigma = Sigma)
  alpha_j            <- mean_satisfaction_score + u[, 1]
  wait_time_effect_j <- mean_wait_time_effect + u[, 2]

  # Make the Linear Predictor 
  mu_satisfaction <- alpha_j[clinic_id]       + 
    wait_time_effect_j[clinic_id] * wait_time +
    beta_emergency * emergency                +
    beta_visits    * visit_num

  # Sample the Satisfaction Scores from Normal Distribution 
  satisfaction_score <- rnorm(n, mean = mu_satisfaction, sd = sigma_satisfaction)

  # Combine in one dataset 
  sim_data <- data.frame(
    satisfaction_score = satisfaction_score,
    wait_time          = wait_time,
    is_emergency       = emergency,
    visits             = visit_num
  )
}
