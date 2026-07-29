
#' Simulate Gym Membership Attendance Data
#'
#' @description
#' ## Scenario
#'
#' A chain of gyms wants to understand what influences how frequently members
#' visit each month.
#'
#' For every member, you observe:
#'
#' * age,
#' * distance from the gym,
#' * whether they work with a personal trainer, and
#' * the total number of gym visits during one month.
#'
#' Members are assumed to be randomly assigned to a trainer, so there is no
#' confounding between trainer assignment and gym attendance.
#'
#' However, gyms are not identical. Different gyms attract different types of
#' members and may have different coaching quality, facilities, or management.
#' Consequently, each gym has its own baseline attendance rate and its own
#' effects of age, travel distance, and trainer assignment.
#'
#' The challenge is to estimate how age, travel distance, and personal trainers
#' influence gym attendance while accounting for variation between gyms.
#'
#' ## Data Generating Model
#'
#' For member \eqn{i} in gym \eqn{j},
#'
#' \deqn{
#' Visits_{ij} \sim Poisson(\lambda_{ij})
#' }
#'
#' where
#'
#' \deqn{
#' \log(\lambda_{ij}) =
#' \alpha_j
#' + \beta_{Age,j} Age_{ij}
#' + \beta_{Distance,j} Distance_{ij}
#' + \beta_{Trainer,j} Trainer_{ij}
#' }
#'
#' Each gym receives its own intercept and regression coefficients:
#'
#' \deqn{
#' \alpha_j
#' \sim
#' Normal(\alpha_0,\sigma_\alpha)
#' }
#'
#' \deqn{
#' \beta_{k,j}
#' \sim
#' Normal(\beta_k,\sigma_k)
#' }
#'
#' where
#' \eqn{k \in \{Age, Distance, Trainer\}}.
#'
#' This simulator therefore generates a varying-intercept,
#' varying-slope hierarchical Poisson regression model.
#' @param gyms Integer. Number of gyms to simulate.
#'
#' @param members_per_gyms Integer. Number of gym members generated for each
#' gym.
#'
#' @param mean_members_age Numeric. Average age (years) of gym members across
#' all gyms.
#'
#' @param prob_has_trainer Numeric. Probability that a member works with a
#' personal trainer.
#'
#' @param mean_log_baseline Numeric. Population-average log expected monthly
#' visit rate before accounting for age, distance, and trainer assignment.
#' For example, \code{log(10)} corresponds to approximately 10 visits per month.
#'
#' @param sigma_gym_visits Numeric. Standard deviation of the gym-specific
#' baseline log visit rates. Larger values produce greater heterogeneity in
#' attendance between gyms.
#'
#' @param mean_age_effect Numeric. Average effect of a one-year increase in age
#' on the log expected number of monthly visits. Negative values imply that
#' older members tend to visit less frequently.
#'
#' @param sigma_age_effects Numeric. Standard deviation of the gym-specific age
#' effects.
#'
#' @param mean_distance_effect Numeric. Average effect of each additional
#' kilometre between a member's home and the gym on the log expected number of
#' monthly visits. Negative values imply that members living farther away visit
#' less frequently.
#'
#' @param sigma_distance_effects Numeric. Standard deviation of the
#' gym-specific distance effects.
#'
#' @param mean_trainer_effect Numeric. Average effect of working with a personal
#' trainer on the log expected number of monthly visits. Positive values imply
#' that members with a trainer attend more frequently.
#'
#' @param sigma_trainer_effects Numeric. Standard deviation of the gym-specific
#' trainer effects.
#'
#' @return A data frame with one row per gym member containing:
#'
#' \describe{
#'   \item{gym_id}{Unique identifier of the gym.}
#'   \item{age}{Age of the gym member (years).}
#'   \item{distance_from_gym}{Distance between the member's home and the gym (km).}
#'   \item{has_trainer}{Indicator for whether the member works with a personal trainer (0 = no, 1 = yes).}
#'   \item{visit}{Number of gym visits during the month.}
#' }
#' @export
simulate_gym_memberships <- function(
  gyms             = 5,
  members_per_gyms = 50,
  mean_members_age = 30,
  prob_has_trainer   = 0.2,
  mean_log_baseline      = log(10),    # ~10 visits/month baseline
  sigma_gym_visits       = 0.20,
  mean_age_effect        = -0.01,      # per year:  e^-0.01 = 0.99x
  sigma_age_effects      = 0.004,
  mean_distance_effect   = -0.08,      # per km:    e^-0.08 = 0.92x
  sigma_distance_effects = 0.02,
  mean_trainer_effect    = 0.25,       # e^0.25 = 1.28x  (28% more visits)
  sigma_trainer_effects  = 0.08

){
  # Simulate Age, Trainers and Distance from the Gym DAG roots
  n       <- gyms * members_per_gyms 
  gyms_id <- rep(seq_len(gyms), each = members_per_gyms) 
  age     <- rnorm(n, mean = mean_members_age, sd = 10)

  # Assuming every gym has a simular training program and aveliability 
  trainer  <- rbinom(n, size = 1, prob = prob_has_trainer) 

  # Distance from the gym is sampled from Uniform distribution in km 
  distance <- runif(n, min = 0.5, max = 5)

  ## Simulation paramters 
  # Baseline Alpha paramter for visits 
  alpha_gym_j <- mean_log_baseline  + sigma_gym_visits * rnorm(n = gyms, mean = 0, sd = 1)
  
  # Age effect on visits per gym 
  age_gym_j <- mean_age_effect + sigma_age_effects * rnorm(n = gyms, mean = 0, sd = 1)  

  # Distance effect on visits per gym 
  distance_gym_j <- mean_distance_effect + sigma_distance_effects * rnorm(n = gyms, mean = 0, sd = 1)

  # Trainer effect on visits per gym 
  trainer_gym_j <- mean_trainer_effect + sigma_trainer_effects * rnorm(n = gyms, mean = 0, sd = 1)

  ## Asseble the Linear Predictor and Sample from Poisson Distribution 
  lambda <- alpha_gym_j[gyms_id]     + 
    age_gym_j[gyms_id]      * age      + 
    distance_gym_j[gyms_id] * distance + 
    trainer_gym_j[gyms_id]  * trainer
  
  visits <- rpois(n = n, lambda = exp(lambda))

  # Combine into one datasets 
  sim_data <- data.frame(
    gym_id = gyms_id,
    age     = age,
    distance_from_gym = distance,
    has_trainer       = trainer,
    visit             = visits
  )
}