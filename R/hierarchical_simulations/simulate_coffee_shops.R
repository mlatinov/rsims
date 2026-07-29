
#' Simulate Coffee Shop Rating Data
#'
#' Generates a hierarchical dataset of coffee shops and baristas.
#' Each barista belongs to a coffee shop, and ratings are generated
#' from a multilevel causal model where shop-level differences,
#' barista experience, and training hours influence customer ratings.
#'
#' The data-generating process is:
#'
#' \deqn{
#' \alpha_j \sim Normal(\mu_{\alpha}, \sigma_{\alpha})
#' }
#'
#' where \eqn{\alpha_j} represents the baseline rating quality of
#' coffee shop \eqn{j}.
#'
#' Barista experience is generated as:
#'
#' \deqn{
#' Experience_i \sim Normal(\mu_E, \sigma_E)
#' }
#'
#' Training hours depend causally on experience:
#'
#' \deqn{
#' Training_i \sim Normal(\mu_T + \beta_T Experience_i,\sigma_T)
#' }
#'
#' Finally, observed ratings are generated as:
#'
#' \deqn{
#' Rating_i \sim Normal(
#' \alpha_j +
#' \beta_E Experience_i +
#' \beta_T Training_i,
#' \sigma_R
#' )
#' }
#'
#' Ratings are restricted to the interval 0,10 .
#'
#' @param num_shops Number of coffee shops.
#' @param baristas_per_shop Number of baristas per coffee shop.
#' @param mean_experience Mean years of barista experience.
#' @param mean_training_hours Baseline mean training hours.
#' @param shop_baseline_score Overall mean coffee shop rating.
#' @param shop_score_sd Standard deviation of shop-level effects.
#' @param experience_training_effect Effect of experience on assigned training hours.
#' @param experience_rating_effect Effect of experience on ratings.
#' @param training_rating_effect Effect of training hours on ratings.
#' @param rating_sd Residual standard deviation of ratings.
#'
#' @return A data.frame containing:
#' \itemize{
#'   \item shop_id: Coffee shop identifier.
#'   \item experience: Barista experience in years.
#'   \item training_hours: Assigned training hours.
#'   \item rating: Simulated customer rating.
#' }
#'
#' @export
simulate_coffee_shops <- function(
  num_shops = 40,
  baristas_per_shop = 4,
  mean_experience     = 2,
  mean_training_hours = 2,
  shops_baseline_score = 5.5,
  shops_score_sd       = 0.5,
  experience_training_effect = 0.1,
  experience_rating_effect   = 0.4,
  training_rating_effect     = 0.2,
  rating_sd = 1
){
  ## Simulate J number of Coffee Shops each with i number of baristas
  n           <- num_shops * baristas_per_shop
  shop_id     <- rep(x = 1:num_shops, each = baristas)
  shops_alpha <- shops_baseline_score + shops_score_sd * rnorm(n = num_shops, mean = 0, sd = 0.3)

  ## Causal Structure and Covariates Experiance Root
  experience <- pmax(rnorm(n, mean = baristas_mean_experiance, sd = 1), 0)     # Barista Experiance in Years 

  # Experience Causes more training hourse to be assignt
  traning_h <- pmax(rnorm(n, mean = baristas_mean_training_hours + experiance_effect_on_training_hours * experience, sd = 0.2))

  # Simulate Ratings from Normal Distribution  
  mu_ratings <- shops_alpha[shop_id] + experience_training_effect * experience + training_rating_effect * traning_h
  ratings    <- pmin(pmax(rnorm(n = n, mean = mu_ratings, sd = ratings_variation),0),10)

  # Return Dataset
  sim_barista <- data.frame(
    shop_id    = shop_id,
    experience = experience,
    traning_hours = traning_h,
    ratings       = ratings
  )
}
