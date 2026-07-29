#' Simulate Exponential Bacterial Colony Growth Using an ODE
#'
#' @description
#' ## Scenario
#'
#' A microbiology experiment follows multiple bacterial cultures during the
#' exponential growth phase, before nutrients become limiting and population
#' growth reaches a ceiling.
#'
#' Each culture is measured daily. Environmental conditions vary between
#' cultures, including nutrient concentration and temperature. These factors
#' influence the bacterial growth rate.
#'
#' The simulator generates hierarchical bacterial growth trajectories where
#' cultures differ in their initial colony size and exponential growth rate.
#'
#' ## Data Generating Model
#'
#' Bacterial growth follows the exponential growth ODE:
#'
#' \deqn{
#' \frac{dC(t)}{dt}=kC(t)
#' }
#'
#' with analytical solution:
#'
#' \deqn{
#' C(t)=C_0\exp(kt)
#' }
#'
#' where \eqn{C_0} is the initial colony size and \eqn{k} is the culture-specific
#' growth rate.
#'
#' Initial colony size is modeled hierarchically:
#'
#' \deqn{
#' \log(C_{0,j})=
#' \alpha_0+\omega_0z_{0,j}
#' }
#'
#' where \eqn{z_{0,j}} represents between-culture variation.
#'
#' The growth rate depends on environmental conditions:
#'
#' \deqn{
#' \log(k_j)=
#' \alpha_k+
#' \beta_T(T_j-\bar{T})
#' +
#' \beta_N(N_j-\bar{N})
#' +
#' \omega_k z_{k,j}
#' }
#'
#' where temperature and nutrient concentration modify the growth rate.
#'
#' Observed bacterial counts include multiplicative measurement error:
#'
#' \deqn{
#' C_{obs}(t)
#' \sim
#' LogNormal(\log(C(t)),\sigma_{obs})
#' }
#'
#' This produces realistic positive bacterial counts with proportional noise.
#'
#' @param num_cultures Integer. Number of independent bacterial cultures.
#'
#' @param num_days Integer. Number of daily measurements per culture.
#'
#' @param mean_nutrient_c Numeric. Population mean nutrient concentration.
#'
#' @param sd_nutrient_c Numeric. Standard deviation of nutrient concentration
#' between cultures.
#'
#' @param mean_temperature Numeric. Population mean culture temperature.
#'
#' @param sd_temperature Numeric. Standard deviation of temperature between
#' cultures.
#'
#' @param alpha_0 Numeric. Population mean log initial colony size.
#'
#' @param omega_0 Numeric. Between-culture standard deviation of log initial
#' colony size.
#'
#' @param alpha_k Numeric. Population mean log growth rate at average
#' temperature and nutrient concentration.
#'
#' @param omega_k Numeric. Between-culture standard deviation of log growth
#' rates.
#'
#' @param beta_temp Numeric. Effect of centered temperature on log growth rate.
#'
#' @param beta_nutrient_c Numeric. Effect of centered nutrient concentration on
#' log growth rate.
#'
#' @param sigma_obs Numeric. Log-scale observation error controlling
#' proportional measurement noise.
#'
#' @return
#' A data frame containing one row per culture-day observation.
#'
#' \describe{
#' \item{culture_id}{Identifier for the bacterial culture.}
#' \item{day}{Measurement day.}
#' \item{nutrient_c}{Nutrient concentration for the culture.}
#' \item{temperature}{Culture temperature.}
#' \item{bacterial_count}{Observed bacterial colony size.}
#' }
#'
#' @examples
#' sim_data <- simulate_ode_bacterial_growth(
#'   num_cultures = 30,
#'   num_days = 20
#' )
#'
#' plot(
#'   bacterial_count ~ day,
#'   data = sim_data,
#'   col = culture_id
#' )
#'
#' @export
simulate_ode_bacterial_growth <- function(
  num_cultures = 30,
  num_days     = 20,
  mean_nutrient_c = 2,   sd_nutrient_c = 0.5,
  mean_temperature = 30, sd_temperature = 2,
  alpha_0 = log(50),     omega_0 = 0.15,   # log C0: intercept + between-culture sd
  alpha_k = log(0.15),  omega_k = 0.10,   # log k:  intercept + between-culture sd
  beta_temp = 0.03,
  beta_nutrient_c = 0.15,
  sigma_obs = 0.08                        # proportional (lognormal) obs error
){
  # Simulate the conditions 
  n   <- num_cultures * num_days
  ids <- make_nested_ids(levels = list(cultures = num_cultures, day = num_days))

  # Simulate the covariates 
  nutrient_c  <- rnorm(num_cultures, mean_nutrient_c, sd_nutrient_c)
  temperature <- rnorm(num_cultures, mean_temperature, sd_temperature)

  # center covariates so alpha_k means "growth rate at average temp & nutrient"
  temp_c <- temperature - mean_temperature
  nutr_c <- nutrient_c  - mean_nutrient_c
  z_0 <- rnorm(num_cultures)
  z_k <- rnorm(num_cultures)

  # Parameters 
  log_C0 <- alpha_0 + omega_0 * z_0
  log_k  <- alpha_k + beta_temp * temp_c + beta_nutrient_c * nutr_c + omega_k * z_k
  C0 <- exp(log_C0)
  k  <- exp(log_k)

  # log C(t) = log C0 + k*t  — do the ODE solution in log space, exponentiate once
  log_mu <- log_C0[ids$cultures_id] + k[ids$cultures_id] * ids$day_id
  bacterial_count <- rlnorm(n, meanlog = log_mu, sdlog = sigma_obs)

  # Return a Dataframe 
  data.frame(
    culture_id  = ids$cultures_id,
    day         = ids$day_id,
    nutrient_c  = nutrient_c[ids$cultures_id],
    temperature = temperature[ids$cultures_id],
    bacterial_count = bacterial_count
  )
}