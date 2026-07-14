#' Simulate Bakery Daily Revenue Data
#'
#' @description
#' ## Scenario
#'
#' A bakery chain wants to understand what drives daily revenue across its
#' stores.
#'
#' For every bakery, daily revenue is recorded together with:
#'
#' * the outdoor temperature (°C),
#' * whether the observation was recorded on a weekend,
#' * the amount of pedestrian foot traffic near the bakery, and
#' * the bakery identifier.
#'
#' At first glance it appears that foot traffic alone explains most of the
#' variation in revenue. However, foot traffic is itself influenced by both
#' the day of the week and the weather. Warm days encourage people to spend
#' more time outside, while weekends change the number of pedestrians passing
#' the bakery.
#'
#' Your task is to estimate the causal effect of foot traffic on bakery revenue.
#' Can you recover the correct effect without adjusting for temperature and
#' weekend status? How do the estimates change after adjustment?
#'
#' This simulator illustrates a classic confounding problem in multiple
#' regression.
#'
#' ## Causal Structure
#'
#' \preformatted{
#'
#' Temperature ───────► Revenue
#'      │
#'      ▼
#' Foot Traffic ─────► Revenue
#'
#' Weekend ──────────► Revenue
#'     │
#'     ▼
#' Foot Traffic
#'
#' }
#'
#' Temperature and weekend status are common causes of both foot traffic
#' and revenue, making them confounders of the relationship between
#' foot traffic and revenue.
#'
#' ## Statistical Task
#'
#' Estimate the effect of pedestrian foot traffic on daily revenue while
#' accounting for temperature, weekend effects, and differences between
#' bakeries. Compare the naïve regression with the correctly adjusted model.
#'
#' ## Data Generating Model
#'
#' Foot traffic is generated as
#'
#' \deqn{
#' FootTraffic_i =
#' \gamma_0
#' +
#' \gamma_1 Weekend_i
#' +
#' \gamma_2 Temperature_i
#' }
#'
#' Daily revenue is then generated as
#'
#' \deqn{
#' Revenue_{ij}
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
#' \beta_T Temperature_{ij}
#' +
#' \beta_F FootTraffic_{ij}
#' +
#' \beta_W Weekend_{ij}
#' }
#'
#' with bakery-specific intercepts
#'
#' \deqn{
#' \alpha_j
#' \sim
#' Normal(\alpha_0,\sigma_\alpha)
#' }
#' @param num_bakeries Integer. Number of bakeries to simulate.
#'
#' @param num_days Integer. Number of daily observations generated for each
#' bakery.
#'
#' @param mean_temperature Numeric. Average outdoor temperature (°C) across
#' all simulated days.
#'
#' @param sd_temperature Numeric. Standard deviation of daily temperatures.
#'
#' @param prob_weekend Numeric. Probability that an observation corresponds
#' to a weekend day.
#'
#' @param mean_foot_traffic Numeric. Average baseline pedestrian foot traffic
#' before accounting for weather and weekends.
#'
#' @param weekend_effect_traffic Numeric. Change in expected foot traffic on
#' weekends relative to weekdays.
#'
#' @param temperature_effect_traffic Numeric. Increase in expected foot traffic
#' for every one-degree increase in outdoor temperature.
#'
#' @param baseline_daily_revenue Numeric. Average baseline daily revenue shared
#' across bakeries before adding bakery-specific variation.
#'
#' @param bakery_revenue_sd Numeric. Standard deviation of bakery-specific
#' baseline revenue.
#'
#' @param beta_temp Numeric. Direct effect of a one-degree increase in
#' temperature on daily revenue, holding all other variables constant.
#'
#' @param beta_traffic Numeric. Direct effect of one additional unit of foot
#' traffic on daily revenue.
#'
#' @param beta_weekend Numeric. Direct effect of weekends on revenue after
#' accounting for foot traffic and temperature.
#'
#' @param daily_revenue_sd Numeric. Residual standard deviation of daily revenue.
#' @export 
simulate_bakeries <- function(
  num_bakeries = 30,
  num_days     = 12,
  mean_temperature = 25,
  sd_temperature   = 8,
  prob_weekend = 0.28,
  mean_foot_trafic = 100,
  weekend_effect_traffic = -40,
  temperature_effect_traffic    = 0.8,
  baseline_daily_revenue = 100,
  bakery_revenue_sd = 20,
  beta_temp   = 0.2,
  beta_trafic = 0.8,
  beta_weekend = 0.4,
  daily_revenue_sd = 10
){

  # Simulate Diffrent Bakaries each having N days of recored daily revenue 
  n <- num_bakeries * num_days
  bakeries_id <- rep(seq_len(num_bakeries), each = num_days)

  # Simulate the Covariates Temperature in C, Foot Trafic in hundereds, is weeked indicator 
  temperature_c <- rnorm(n, mean = mean_temperature, sd = sd_temperature)
  is_weekend     <- rbinom(n, size = 1, prob = prob_weekend)

  # The foot Trafic is partially caused by bolth if is a weekend and the outside temperature
  foot_traffic <- mean_foot_trafic + weekend_effect_traffic * is_weekend + temperature_effect_traffic * temperature_c 

  # Make a Linear Predictor 
  revenue_baseline <- baseline_daily_revenue + bakeries_revenue_deviation * rnorm(num_bakeries, mean = 0, sd = 1)
  daily_revenue_mu <- revenue_baseline[bakeries_id] 
    + beta_temp * temperature_c
    + beta_trafic  * foot_traffic 
    + beta_weekend  * is_weekend

  # Sample from normal Distribution 
  daily_revenue <- rnorm(n, mean = daily_revenue_mu, sd = daily_revenue_sd)

  # Combine in one dataset 
  sim_data <- data.frame(
    daily_revenue = daily_revenue,
    foot_traffic   = foot_traffic,
    is_weekend    = is_weekend,
    temperature_c = temperature_c,
    bakery_id    = bakeries_id
  )
}
