#' Simulate Restaurant Ratings with Crossed Restaurant and Cuisine Effects
#'
#' @description
#' ## Scenario
#'
#' A food review company wants to understand which factors influence customer
#' ratings of restaurant dishes.
#'
#' Every restaurant offers every cuisine in the study, creating a crossed
#' experimental design where cuisines are evaluated across multiple
#' restaurants.
#'
#' For every restaurant-cuisine combination, you observe:
#'
#' * whether the dish is the restaurant's signature (special) dish,
#' * the dish price,
#' * the preparation time,
#' * the customer rating, and
#' * the restaurant and cuisine identifiers.
#'
#' Signature dishes generally cost more and require longer preparation times,
#' but customers also tend to rate them more favorably.
#'
#' Restaurants differ in overall service quality, atmosphere, and staff,
#' whereas cuisines differ in their intrinsic popularity and appeal. Both
#' sources of variation contribute independently to customer ratings.
#'
#' The challenge is to estimate the effects of price, waiting time, and
#' signature dishes while accounting for crossed restaurant and cuisine
#' effects.
#'
#' ## Causal Structure
#'
#' \preformatted{
#'
#'                  Signature Dish
#'                    │      │
#'                    │      ├──────────────► Price ───────────────┐
#'                    │      │                                     │
#'                    │      └──────────────► Waiting Time ────────┤
#'                    │                                            ▼
#'                    └──────────────────────────────────────────► Rating
#'
#' Restaurant ───────────────────────────────────────────────────► Rating
#'
#' Cuisine ──────────────────────────────────────────────────────► Rating
#'
#' }
#'
#' Restaurant and cuisine effects are crossed rather than nested.
#'
#' ## Statistical Task
#'
#' Estimate the effects of price, preparation time, and signature dishes on
#' customer ratings while accounting for both restaurant-specific and
#' cuisine-specific variation. Compare a standard linear regression with a
#' crossed-effects multilevel model including varying intercepts for both
#' restaurants and cuisines.
#'
#' ## Data Generating Model
#'
#' Dish prices are generated as
#'
#' \deqn{
#' Price_i
#' \sim
#' Normal(\mu_i,\sigma_{Price})
#' }
#'
#' where
#'
#' \deqn{
#' \mu_i
#' =
#' \alpha_P
#' +
#' \beta_{Special}Special_i
#' }
#'
#' Waiting times are generated as
#'
#' \deqn{
#' Wait_i
#' \sim
#' Poisson(\lambda_i)
#' }
#'
#' where
#'
#' \deqn{
#' \lambda_i
#' =
#' \alpha_W
#' +
#' \beta_{Special}Special_i
#' }
#'
#' Customer ratings are generated as
#'
#' \deqn{
#' Rating_{ij}
#' \sim
#' Normal(\mu_{ij},\sigma)
#' }
#'
#' where
#'
#' \deqn{
#' \mu_{ij}
#' =
#' \alpha
#' +
#' \alpha_j^{Restaurant}
#' +
#' \alpha_k^{Cuisine}
#' +
#' \beta_{Price}Price_{ij}
#' +
#' \beta_{Wait}Wait_{ij}
#' +
#' \beta_{Special}Special_{ij}
#' }
#'
#' Restaurant and cuisine intercepts are independently generated as
#'
#' \deqn{
#' \alpha_j^{Restaurant}
#' \sim
#' Normal(0,\sigma_{Restaurant})
#' }
#'
#' \deqn{
#' \alpha_k^{Cuisine}
#' \sim
#' Normal(0,\sigma_{Cuisine})
#' }
#'
#' @param num_restaurants Integer. Number of restaurants.
#'
#' @param num_cuisines Integer. Number of cuisines.
#'
#' @param baseline_dish_price Numeric. Average baseline dish price.
#'
#' @param p_beta_special Numeric. Effect of a signature dish on expected price.
#'
#' @param sd_price Numeric. Residual standard deviation of dish prices.
#'
#' @param baseline_wait Numeric. Baseline preparation time.
#'
#' @param w_beta_special Numeric. Effect of a signature dish on expected
#' preparation time.
#'
#' @param baseline_rating_score Numeric. Population-average customer rating.
#'
#' @param sd_restaurants_baseline_score Numeric. Standard deviation of
#' restaurant-specific intercepts.
#'
#' @param sd_cuisines_baseline_score Numeric. Standard deviation of
#' cuisine-specific intercepts.
#'
#' @param beta_price Numeric. Effect of dish price on customer ratings.
#'
#' @param beta_wait Numeric. Effect of preparation time on customer ratings.
#'
#' @param beta_special Numeric. Direct effect of signature dishes on customer
#' ratings.
#'
#' @param sd_rating_score Numeric. Residual standard deviation of ratings.
#'
#' @return A data frame with one row for every restaurant-cuisine combination.
#'
#' \describe{
#' \item{restaurant_id}{Restaurant identifier.}
#' \item{cuisine_id}{Cuisine identifier.}
#' \item{is_special}{Indicator for whether the dish is the restaurant's signature dish.}
#' \item{price}{Dish price.}
#' \item{wait_time}{Preparation time.}
#' \item{rating}{Customer rating on a 0–100 scale.}
#' }
#'
#' @export
simulate_restaurants_x_cuisines <- function(
  num_restaurants = 8,
  num_cuisines = 25,
  baseline_dish_price = 18,
  p_beta_special = 6,
  sd_price = 2,
  baseline_wait = 12,
  w_beta_special = 4,
  baseline_rating_score = 75,
  sd_restaurants_baseline_score = 4,
  sd_cuisines_baseline_score = 5,
  beta_price = -0.35,
  beta_wait = -0.50,
  beta_special = 8,
  sd_rating_score = 6
){

  # Simuate the Crossed Desing between estaurants and cuisines 
  n <- num_restaurants * num_cuisines
  restaurants_id <- rep(seq_len(num_restaurants), each  = num_cuisines)
  cuisines_id    <- rep(seq_len(num_cuisines),    times = num_restaurants)

  # Simulate the Root Variable Is the dish the restourant special 
  is_special <- rbinom(n, size = 1, prob = 0.2)

  # Simulate the price of the dish The special dish effect the price tag
  mu_price <- baseline_dish_price + p_beta_special * is_special  
  price    <- pmax(rnorm(n, mean = mu_price, sd = sd_price), 5) 

  # Simulate the Preparation wait time as dependent on the is the dish special 
  mu_wait   <- baseline_wait + w_beta_special * is_special
  wait_time <- rpois(n, lambda = mu_wait) 

  ## Simulate Rating Scores 1-100 
  # Rating Scores Paramters 
  alpha_restorant <- sd_restaurants_baseline_score * rnorm(num_restaurants, mean = 0, sd = 1) 
  alpha_cuisines  <- sd_cuisines_baseline_score * rnorm(num_cuisines, mean = 0, sd = 1) 

  mu_rating <- (
    baseline_rating_score 
      + alpha_restorant[restaurants_id] 
      + alpha_cuisines[cuisines_id]
      + beta_price   * price
      + beta_wait    * wait_time
      + beta_special * is_special
    )
  
  # Sample the rating scores from Normal Distribution 
  ratings <- pmin(pmax(rnorm(n, mean = mu_rating, sd = sd_rating_score),0), 100)

  # Combine in one simulated dataset 
  sim_data <- data.frame(
    restaurant_id = restaurants_id,
    cuisines_id   = cuisines_id,
    price         = price,
    wait_time     = wait_time,
    is_special    = is_special,
    rating        = ratings 
  )
}
