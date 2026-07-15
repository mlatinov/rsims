#' Simulate Hierarchical Hospital Recovery Data
#'
#' @description
#' ## Scenario
#'
#' A healthcare network wants to evaluate whether a standardized treatment dose
#' improves patient recovery across multiple hospitals.
#'
#' Patients are admitted to hospital wards, and wards are nested within
#' hospitals.
#'
#' For every patient, you observe:
#'
#' * standardized treatment dose,
#' * age,
#' * number of comorbidities,
#' * whether the patient recovered, and
#' * the hospital and ward identifiers.
#'
#' Older patients tend to have more comorbidities, which reduce the probability
#' of recovery. Hospitals differ in their baseline recovery rates because of
#' differences in equipment, staff expertise, and patient management.
#' Furthermore, hospitals differ in how effective the treatment dose is,
#' producing hospital-specific treatment effects.
#'
#' Wards within hospitals also differ in their baseline recovery rates because
#' of differences in staffing, patient mix, or organization.
#'
#' The challenge is to estimate the effect of treatment dose while accounting
#' for patient-level covariates and the hierarchical hospital structure.
#'
#' ## Causal Structure
#'
#' \preformatted{
#'
#' Age ───────────────► Comorbidities ─────────────► Recovery
#'  │                                               ▲
#'  └───────────────────────────────────────────────┘
#'
#' Treatment Dose ─────────────────────────────────► Recovery
#'
#' Hospital ─────────► Ward ───────────────────────► Recovery
#'        │
#'        └────────────────────────────────────────► Dose Effect
#'
#' }
#'
#' Hospitals have varying baseline recovery rates and varying treatment
#' effects, while wards contribute additional nested variation.
#'
#' ## Statistical Task
#'
#' Estimate the effect of treatment dose on patient recovery while accounting
#' for age, comorbidities, hospital-level variation, and ward-level variation.
#' Compare a standard logistic regression with a hierarchical logistic model
#' including varying intercepts for hospitals and wards and varying treatment
#' effects across hospitals.
#'
#' ## Data Generating Model
#'
#' Comorbidities are generated according to
#'
#' \deqn{
#' Comorbidities_i
#' \sim
#' Poisson(\lambda_i)
#' }
#'
#' where
#'
#' \deqn{
#' \log(\lambda_i)
#' =
#' \alpha_C
#' +
#' \beta_{Age}Age_i
#' }
#'
#' Patient recovery is generated as
#'
#' \deqn{
#' Recovery_i
#' \sim
#' Bernoulli(\pi_i)
#' }
#'
#' where
#'
#' \deqn{
#' \mathrm{logit}(\pi_i)
#' =
#' \alpha_h
#' +
#' \alpha_w
#' +
#' \beta_{Dose,h}Dose_i
#' +
#' \beta_{Comorbidity}Comorbidities_i
#' +
#' \beta_{Age}Age_i
#' }
#'
#' Hospital-specific intercepts and treatment effects are jointly generated as
#'
#' \deqn{
#' (\alpha_h,\beta_{Dose,h})
#' \sim
#' MVN(\mathbf{0},\Sigma)
#' }
#'
#' while ward intercepts follow
#'
#' \deqn{
#' \alpha_w
#' \sim
#' Normal(0,\sigma_{Ward})
#' }
#'
#' @param num_hospitals Integer. Number of hospitals.
#' @param wards_per_hospital Integer. Number of wards within each hospital.
#' @param num_patients_per_ward Integer. Number of patients treated in each ward.
#' @param baseline_comorbidity_rate Numeric. Baseline log-rate of comorbidities.
#' @param c_beta_age Numeric. Effect of age on the expected number of comorbidities.
#' @param rho Numeric. Correlation between hospital baseline recovery and treatment effect.
#' @param sd_hospital_level_alpha Numeric. Standard deviation of hospital-specific intercepts.
#' @param sd_dose_level_hospital_beta Numeric. Standard deviation of hospital-specific dose effects.
#' @param baseline_recovery_hospital Numeric. Population-average hospital intercept.
#' @param baseline_dose_effect Numeric. Average treatment-dose effect.
#' @param baseline_wards_recovery Numeric. Population-average ward intercept.
#' @param sd_wards_recovery Numeric. Standard deviation of ward intercepts.
#' @param beta_comorbidities Numeric. Effect of each additional comorbidity on recovery.
#' @param beta_age Numeric. Direct effect of patient age on recovery.
#'
#' @return A data frame with one row per patient.
#'
#' \describe{
#' \item{hospital_id}{Hospital identifier.}
#' \item{ward_id}{Ward identifier nested within hospitals.}
#' \item{std_dose}{Standardized treatment dose.}
#' \item{age}{Patient age in years.}
#' \item{comorbidities}{Number of diagnosed comorbidities.}
#' \item{recovery}{Recovery outcome (0 = no, 1 = yes).}
#' }
#'
#' @export
simulate_hospitals_wards <- function(
  num_hospitals = 20,
  wards_per_hospital = 4,
  num_patients_per_ward = 20,

  # Comorbidity model
  baseline_comorbidity_rate = -1.8,
  c_beta_age = 0.03,

  # Correlation between hospital intercept and treatment effect
  rho = 0.55,

  # Hospital-level variation
  sd_hospital_level_alpha = 0.4,
  sd_dose_level_hospital_beta = 0.25,

  # Population-average effects
  baseline_recovery_hospital = -0.4,
  baseline_dose_effect = 0.8,

  # Ward-level variation
  baseline_wards_recovery = 0,
  sd_wards_recovery = 0.35,

  # Patient-level effects
  beta_comorbidities = -0.35,
  beta_age = -0.015,

  # Residual variation
  dose_mean = 0,
  dose_sd = 1
){

  # Simulate J hospitals each containing K wards with N patients per ward
  n <- num_hospitals * wards_per_hospital * num_patients_per_ward
  hospital_id <- rep(rep(seq_len(num_hospitals), each = wards_per_hospital),each = num_patients_per_ward)
  ward_id <- rep(seq_len(num_hospitals * wards_per_hospital),each = num_patients_per_ward)

  # Simulate standardized treatment dose
  dose <- rnorm(n, mean = dose_mean, sd = dose_sd)

  # Simulate patient age
  age <- rnorm(n, mean = 60, sd = 10)

  # Simulate comorbidities dependent on age
  lambda_comorbidities <- exp(baseline_comorbidity_rate + c_beta_age * age)
  comorbidities <- rpois(n, lambda = lambda_comorbidities)

  # Correlated hospital-specific intercepts and dose effects
  cor_matrix <- matrix(c(1, rho, rho, 1),nrow = 2, byrow = TRUE)
  sigmas <- diag(c(sd_hospital_level_alpha, sd_dose_level_hospital_beta))
  Sigma <- sigmas %*% cor_matrix %*% sigmas

  # Sample hospital-level effects
  u <- MASS::mvrnorm(n = num_hospitals, mu = c(0, 0), Sigma = Sigma)
  alpha_h   <- baseline_recovery_hospital + u[, 1]
  beta_dose <- baseline_dose_effect + u[, 2]

  # Sample ward-level intercepts
  alpha_w <- baseline_wards_recovery + sd_wards_recovery * rnorm(num_hospitals * wards_per_hospital, mean = 0, sd = 1)

  # Construct linear predictor
  logit_pi <-
    alpha_h[hospital_id] +
    alpha_w[ward_id] +
    beta_dose[hospital_id] * dose +
    beta_comorbidities * comorbidities +
    beta_age * age

  # Convert to probabilities
  recovery_prob <- plogis(logit_pi)

  # Sample recovery outcome
  recovery <- rbinom(n, size = 1, prob = recovery_prob)

  # Return simulated data
  sim_data <- data.frame(
    hospital_id = hospital_id,
    ward_id = ward_id,
    std_dose = dose,
    age = age,
    comorbidities = comorbidities,
    recovery = recovery
  )
  sim_data
}