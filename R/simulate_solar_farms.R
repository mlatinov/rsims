#' Simulate Solar Farm Energy Production
#'
#' @description
#' ## Scenario
#'
#' An energy company operates multiple solar farms and wants to understand
#' which environmental conditions influence daily electricity production.
#'
#' For every farm and day, you observe:
#'
#' * hours of sunlight,
#' * average solar panel temperature,
#' * dust accumulation on the panels,
#' * daily electricity production (kWh), and
#' * the farm identifier.
#'
#' Longer sunlight exposure generally increases energy production but also
#' raises panel temperature. Hotter panels become less efficient, reducing
#' electricity generation. Dust accumulation blocks incoming sunlight and
#' further decreases production.
#'
#' Solar farms differ in their baseline productivity, their sensitivity to
#' sunlight, and how strongly panel temperature affects efficiency. These
#' farm-specific effects are allowed to be correlated.
#'
#' The challenge is to estimate the direct effects of sunlight, panel
#' temperature, and dust accumulation while accounting for between-farm
#' variation and correlated varying effects.
#'
#' ## Causal Structure
#'
#' \preformatted{
#'
#' Sun Hours ─────────────► Panel Temperature ─────────────► Daily kWh
#'      │                           │
#'      └───────────────────────────┘
#'
#' Dust Index ────────────────────────────────────────────► Daily kWh
#'
#' Farm ──────────────────────────────────────────────────► Daily kWh
#'
#' }
#'
#' Panel temperature acts as a mediator between sunlight exposure and
#' electricity production.
#'
#' ## Statistical Task
#'
#' Estimate the effects of sunlight, panel temperature, and dust accumulation
#' on daily electricity generation while accounting for correlated
#' farm-specific intercepts and slopes. Compare a standard regression with a
#' multilevel model using correlated varying effects.
#'
#' ## Data Generating Model
#'
#' Panel temperature is generated as
#'
#' \deqn{
#' Temperature_i
#' \sim
#' Normal(\mu_i,\sigma_{Temp})
#' }
#'
#' where
#'
#' \deqn{
#' \mu_i
#' =
#' \alpha_{Temp}
#' +
#' \beta_{Sun}SunHours_i
#' }
#'
#' Daily electricity production is generated as
#'
#' \deqn{
#' kWh_{ij}
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
#' \beta_jSunHours_{ij}
#' +
#' \gamma_jPanelTemperature_{ij}
#' +
#' \delta DustIndex_{ij}
#' }
#'
#' The farm-specific effects are jointly distributed as
#'
#' \deqn{
#' \left(
#' \alpha_j,
#' \beta_j,
#' \gamma_j
#' \right)
#' \sim
#' MVN(\mathbf{0},\Sigma)
#' }
#'
#' allowing the intercept and slopes to be correlated across farms.
#'
#' @param num_farms Integer. Number of solar farms.
#'
#' @param days_per_farm Integer. Number of simulated days for each farm.
#'
#' @param sun_hours_received_range Numeric vector giving the minimum and
#' maximum daily sunlight hours.
#'
#' @param baseline_panel_temp Numeric. Baseline panel temperature (°C).
#'
#' @param panel_temp_beta_sun_hours Numeric. Effect of sunlight hours on panel
#' temperature.
#'
#' @param p_sd_panel_temp Numeric. Residual standard deviation of panel
#' temperature.
#'
#' @param baseline_kw Numeric. Population-average daily electricity production.
#'
#' @param sun_hours_effect_kw Numeric. Average effect of one additional hour of
#' sunlight on daily production.
#'
#' @param panel_temp_effect_kw Numeric. Average effect of panel temperature on
#' electricity production.
#'
#' @param dust_effect_kw Numeric. Effect of dust accumulation on electricity
#' production.
#'
#' @param sd_kw Numeric. Residual standard deviation of daily electricity
#' production.
#'
#' @param sd_baseline_kw Numeric. Standard deviation of farm-specific baseline
#' production.
#'
#' @param sd_sun_hours Numeric. Standard deviation of the farm-specific
#' sunlight effects.
#'
#' @param sd_panel_temp Numeric. Standard deviation of the farm-specific panel
#' temperature effects.
#'
#' @param cor_matrix Correlation matrix governing the farm-specific intercept
#' and slope effects.
#'
#' @return A data frame with one row per farm-day observation.
#'
#' \describe{
#' \item{farm_id}{Solar farm identifier.}
#' \item{sun_hours}{Daily sunlight exposure (hours).}
#' \item{panel_temp}{Average panel temperature (°C).}
#' \item{dust_index}{Dust accumulation index between 0 and 1.}
#' \item{kw}{Daily electricity production (kWh).}
#' }
#'
#' @export
simulate_solar_farms <- function(
  num_farms = 25,
  days_per_farm = 120,
  sun_hours_received_range = c(4, 11),
  baseline_panel_temp = 18,
  panel_temp_beta_sun_hours = 2.5,
  p_sd_panel_temp = 2.5,
  baseline_kw = 55,
  sun_hours_effect_kw = 6,
  panel_temp_effect_kw = -0.75,
  dust_effect_kw = -12,
  sd_kw = 4,
  sd_baseline_kw = 6,
  sd_sun_hours = 0.35,
  sd_panel_temp = 0.08,
  cor_matrix = matrix(c(1, -0.4, 0.3, -0.4, 1, -0.2, 0.3, -0.2, 1),nrow = 3, ncol = 3)
){
  # Simulate the N Solar Farms each runed for J number of days
  n       <- num_farms * days_per_farm
  farm_id <- rep(seq_len(num_farms), each = days_per_farm) 

  # Simulate Root Dust Index as value between 0 and 1 
  dust_index <- runif(n, min = 0, max = 1)

  # Simulate Sun Hours reserved 
  sun_hours  <- runif(n, min = sun_hours_received_range[1], max = sun_hours_received_range[2]) 

  # Simulate Panel Temperature dependent on the sun hours exposure 
  mu_panel_temp <- baseline_panel_temp + panel_temp_beta_sun_hours * sun_hours
  panel_temp    <- rnorm(n, mean = mu_panel_temp, sd = p_sd_panel_temp)

  # Simulate Kw with correlated intercept and panel temp and sun_h correlation 
  sigmas <- diag(c(sd_baseline_kw, sd_sun_hours, sd_panel_temp))
  Sigma  <- sigmas %*% cor_matrix %*% sigmas

  # Sample the paramters from MVN distribution 
  u <- MASS::mvrnorm(n = num_farms, mu = c(0, 0, 0), Sigma = Sigma)

  # Compose the effects 
  alpha_j <- baseline_kw + u[,1]
  beta_j  <- sun_hours_effect_kw + u[,2]
  gamma_j <- panel_temp_effect_kw + u[,3]

  # Linear Predictor 
  mu_kw <- (
    alpha_j[farm_id] 
      + beta_j[farm_id]  * sun_hours
      + gamma_j[farm_id] * panel_temp 
      + dust_effect_kw  * dust_index
  )

  # Sample from Kw per day from Normal Distribution 
  kw <- rnorm(n, mean = mu_kw, sd = sd_kw)

  # Combine in simulated dataframe 
  sim_data <- data.frame(
    farm_id    = farm_id,
    sun_hours  = sun_hours,
    panel_temp = panel_temp,
    dust_index = dust_index,
    kw = kw
  )
}
