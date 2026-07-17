#' Simulate Plant Biomass Response to Fertiliser Application
#'
#' @description
#' Simulates a hierarchical nonlinear relationship between fertiliser quantity
#' and plant biomass across multiple soil types.
#'
#' Each soil type is treated as a separate group with its own baseline biomass
#' level and its own nonlinear fertiliser response curve. The fertiliser effect
#' is represented using a hierarchical spline model, allowing soil types to have
#' similar but not identical response functions.
#'
#' The data generating process is
#'
#' \deqn{
#' Biomass_i \sim Normal(\mu_i,\sigma)
#' }
#'
#' with
#'
#' \deqn{
#' \mu_i =
#' \alpha_j + f_j(Fertiliser_i)
#' }
#'
#' where \eqn{j} denotes the soil type.
#'
#' Soil-specific baseline biomass values are generated as
#'
#' \deqn{
#' \alpha_j =
#' \alpha_0 + u_j,
#' \qquad
#' u_j \sim Normal(0,\sigma_\alpha)
#' }
#'
#' The nonlinear fertiliser response is generated using a hierarchical spline:
#'
#' \deqn{
#' f_j(x_i)=B(x_i)w_j
#' }
#'
#' where the soil-specific spline weights follow
#'
#' \deqn{
#' w_j=\bar w+\epsilon_j,
#' \qquad
#' \epsilon_j\sim Normal(0,\sigma_{group})
#' }
#'
#' This creates partially pooled fertiliser response curves where soil types
#' share information while retaining their own nonlinear patterns.
#'
#' @param num_soil_types Number of soil types simulated.
#'
#' @param num_of_plots_per_soil_type Number of experimental plots generated for
#' each soil type.
#'
#' @param baseline_soil_type_biomass Population average baseline biomass before
#' accounting for soil-specific variation.
#'
#' @param soil_type_biomass_sd Standard deviation controlling variation in
#' baseline biomass between soil types.
#'
#' @param fertiliser_weights Population-level spline weights controlling the
#' average nonlinear fertiliser response curve.
#'
#' @param plant_biomass_sd Residual standard deviation of plant biomass around
#' the expected value.
#'
#' @return
#' A data frame containing simulated plant growth measurements:
#'
#' \describe{
#'   \item{soil_type_id}{Identifier of the soil type.}
#'   \item{fertiliser}{Amount of fertiliser applied.}
#'   \item{plant_biomass}{Observed plant biomass.}
#' }
#'
#' @details
#' Fertiliser response is generated using radial basis functions. The basis
#' functions are fixed transformations of fertiliser quantity, while the spline
#' coefficients determine the nonlinear response shape.
#'
#' The supplied `fertiliser_weights` represent the population-average response.
#' Soil-specific deviations are generated through the hierarchical spline helper,
#' producing different response curves for different soil types.
#'
#' @examples
#' set.seed(123)
#'
#' data <- simulate_plant_growth_vs_fertiliser()
#'
#' plot(
#'   data$fertiliser,
#'   data$plant_biomass,
#'   xlab = "Fertiliser quantity",
#'   ylab = "Plant biomass"
#' )
#'
#' @export
simulate_plant_growth_vs_fertiliser <- function(
  num_soil_types = 6,
  num_of_plots_per_soil_type = 100,
  baseline_soil_type_biomass = 3,
  soil_type_biomass_sd       = 1,
  fertiliser_weights = c(2, 9, 13, 8, 1),
  plant_biomass_sd   = 1
){

  # Simulate J number of soil types with N number of plots 
  idx <- make_nested_ids(levels = c(soil_type = num_soil_types, plots = num_of_plots_per_soil_type))

  # Simulate the fertiliser quantity 
  fertiliser <- runif(n = nrow(idx), min = 0, max = 100)

  # Paramaters 
  alpha_j  <- baseline_soil_type_biomass + soil_type_biomass_sd * rnorm(n = num_soil_types, mean = 0, sd = 1)
  B        <- make_radial_basis(x = fertiliser, knots = 5, length_scale = 22)
  f_fertiliser_j <- simulate_hierarchical_spline_effect(
    B = B,
    group       = idx$soil_type_id,
    sigma_group = 2, 
    weights     = fertiliser_weights
  )

  # Make a Linear Predictor and sample from Normal Distribution 
  mu_i    <- (alpha_j[idx$soil_type_id] + f_fertiliser_j$fitted)
  biomass <- rnorm(n = nrow(idx), mean = mu_i, sd = plant_biomass_sd) 

  # Combine in one dataset 
  sim_data <- data.frame(
    soil_type_id  = idx$soil_type_id,
    fertiliser    = fertiliser,
    plant_biomass = biomass 
  )
}  
