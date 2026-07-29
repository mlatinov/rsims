#' Simulate Bacterial Colony Growth Using a Geometric Brownian Motion SDE
#'
#' @description
#' ## Scenario
#'
#' A microbiology laboratory monitors multiple bacterial cultures during the
#' exponential growth phase, before nutrient depletion or carrying capacity
#' effects become relevant.
#'
#' Each culture is measured daily, with environmental conditions recorded:
#'
#' * nutrient concentration,
#' * incubation temperature,
#' * bacterial colony size.
#'
#' Different cultures have different initial colony sizes and intrinsic growth
#' rates. Temperature and nutrient availability influence the expected growth
#' rate, while stochastic fluctuations represent biological variability in the
#' growth process.
#'
#' This simulator extends the deterministic exponential growth ODE into a
#' stochastic differential equation model using geometric Brownian motion (GBM).
#' The process remains positive while allowing random deviations around the
#' expected exponential trajectory.
#'
#' ## Deterministic Growth Model (ODE)
#'
#' Without stochastic variation, bacterial growth follows:
#'
#' \deqn{
#' \frac{dC(t)}{dt}=kC(t)
#' }
#'
#' with solution:
#'
#' \deqn{
#' C(t)=C_0\exp(kt)
#' }
#'
#' where \eqn{C_0} is the initial colony size and \eqn{k} is the culture-specific
#' growth rate.
#'
#' ## Stochastic Growth Model (SDE)
#'
#' The stochastic extension is:
#'
#' \deqn{
#' dC(t)=kC(t)dt+\sigma C(t)dW(t)
#' }
#'
#' where:
#'
#' * \eqn{W(t)} is a Brownian motion process,
#' * \eqn{\sigma} controls biological growth volatility,
#' * the stochastic term introduces multiplicative uncertainty.
#'
#' The simulator uses the exact discrete-time GBM transition:
#'
#' \deqn{
#' C_{t+\Delta t}
#' =
#' C_t
#' \exp
#' \left[
#' (k-\frac{\sigma^2}{2})\Delta t
#' +
#' \sigma\sqrt{\Delta t}Z
#' \right]
#' }
#'
#' where \eqn{Z \sim Normal(0,1)}.
#'
#' ## Hierarchical Growth Parameters
#'
#' Initial colony size varies between cultures:
#'
#' \deqn{
#' \log(C_{0,j})
#' =
#' \alpha_0+\omega_0 z_{0,j}
#' }
#'
#' Growth rates vary between cultures and depend on environmental covariates:
#'
#' \deqn{
#' \log(k_j)
#' =
#' \alpha_k
#' +
#' \beta_{temp}(Temp_j-\bar{Temp})
#' +
#' \beta_{nutrient}(Nutrient_j-\bar{Nutrient})
#' +
#' \omega_k z_{k,j}
#' }
#'
#' where \eqn{z_{0,j}} and \eqn{z_{k,j}} represent culture-specific random
#' effects.
#'
#' @param num_cultures Integer. Number of bacterial cultures simulated.
#'
#' @param num_days Integer. Number of daily measurements collected for each
#' culture.
#'
#' @param mean_nutrient_c Mean nutrient concentration across cultures.
#'
#' @param sd_nutrient_c Standard deviation of nutrient concentration.
#'
#' @param mean_temperature Mean incubation temperature.
#'
#' @param sd_temperature Standard deviation of incubation temperature.
#'
#' @param alpha_0 Population mean log initial colony size.
#'
#' @param omega_0 Standard deviation of culture-specific variation in initial
#' colony size.
#'
#' @param alpha_k Population mean log growth rate.
#'
#' @param omega_k Standard deviation of culture-specific variation in growth
#' rate.
#'
#' @param beta_temp Effect of temperature deviation on bacterial growth rate.
#'
#' @param beta_nutrient_c Effect of nutrient concentration deviation on
#' bacterial growth rate.
#'
#' @param sigma_gbm Diffusion volatility parameter controlling stochastic
#' fluctuations in bacterial growth. This represents process uncertainty and is
#' not measurement error.
#'
#' @param dt Time step size used for the GBM transition.
#'
#' @return
#' A data frame containing one row per culture-day observation.
#'
#' \describe{
#' \item{culture_id}{Identifier of the bacterial culture.}
#' \item{day}{Measurement day.}
#' \item{nutrient_c}{Nutrient concentration associated with the culture.}
#' \item{temperature}{Incubation temperature associated with the culture.}
#' \item{bacterial_count}{Simulated bacterial colony size generated from the GBM
#' process.}
#' }
#'
#' @examples
#' sim_data <- simulate_sde_gbm_bacterial_growth()
#'
#' plot(
#'   sim_data$day,
#'   sim_data$bacterial_count,
#'   type = "l"
#' )
#'
#' @export
simulate_sde_gbm_bacterial_growth <- function(
  num_cultures = 30,
  num_days     = 20,
  mean_nutrient_c = 2,   sd_nutrient_c = 0.5,
  mean_temperature = 30, sd_temperature = 2,
  alpha_0 = log(50),     omega_0 = 0.15,
  alpha_k = log(0.15),   omega_k = 0.10,
  beta_temp = 0.03,
  beta_nutrient_c = 0.15,
  sigma_gbm = 0.08,     # process (diffusion) volatility -- NOT observation noise
  dt = 1
){
  # Simulate Covariates
  nutrient_c  <- rnorm(num_cultures, mean_nutrient_c, sd_nutrient_c)
  temperature <- rnorm(num_cultures, mean_temperature, sd_temperature)

  # Parameters
  temp_c <- temperature - mean_temperature
  nutr_c <- nutrient_c  - mean_nutrient_c
  z_0 <- rnorm(num_cultures)
  z_k <- rnorm(num_cultures)
  log_C0 <- alpha_0 + omega_0 * z_0
  log_k  <- alpha_k + beta_temp * temp_c + beta_nutrient_c * nutr_c + omega_k * z_k
  C0 <- exp(log_C0)
  k  <- exp(log_k)              # length num_cultures -- one rate per culture, reused every step

  # Simulate Geometric Brownian motion with Ito correction 
  C <- matrix(NA_real_, nrow = num_days, ncol = num_cultures)
  C[1, ] <- C0
  for (t in 2:num_days) {
    drift_t <- (k - sigma_gbm^2 / 2) * dt                     
    C[t, ]  <- rlnorm(num_cultures, meanlog = log(C[t-1, ]) + drift_t, sdlog   = sigma_gbm * sqrt(dt))        
  }

  # only flatten to long format AFTER the recursion is done
  culture_id <- rep(1:num_cultures, each = num_days)
  day        <- rep(1:num_days, times = num_cultures)

  data.frame(
    culture_id  = culture_id,
    day         = day,
    nutrient_c  = nutrient_c[culture_id],
    temperature = temperature[culture_id],
    bacterial_count = as.vector(t(C))
  )
}