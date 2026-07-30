#' Simulate Contaminant Diffusion Through Soil with Site-Level Diffusivity
#'
#' @description
#' Simulates a hierarchical contaminant diffusion process across multiple soil
#' sites using a one-dimensional diffusion partial differential equation (PDE).
#'
#' Each site has its own soil porosity value, which influences the site-specific
#' diffusion coefficient:
#'
#' \deqn{
#' \log(D_j) =
#' \alpha_D +
#' \beta_D \times porosity_j +
#' \omega_j
#' }
#'
#' where:
#'
#' * `D_j` is the diffusivity at site `j`,
#' * `porosity_j` is a site-level soil property,
#' * `\omega_j` represents unexplained site-level variation.
#'
#' The contaminant concentration field evolves according to a pure diffusion
#' equation:
#'
#' \deqn{
#' \frac{\partial u_j(x,t)}{\partial t}
#' =
#' D_j
#' \frac{\partial^2 u_j(x,t)}{\partial x^2}
#' }
#'
#' where `u_j(x,t)` represents contaminant concentration at depth `x` and time
#' `t` for site `j`.
#'
#' The simulation separates:
#'
#' * site-level covariates,
#' * PDE parameters,
#' * spatial discretization,
#' * temporal integration,
#' * numerical stability control,
#' * observation sampling.
#'
#' The full concentration field is simulated internally, but measurements are
#' only collected at selected sensor depths and times, mimicking a real soil
#' monitoring experiment.
#'
#' The observation layer follows:
#'
#' \deqn{
#' y_{jkl} \sim LogNormal(
#' \log(u_j(x_l,t_k)),
#' \sigma_{obs}
#' )
#' }
#'
#' where observations are noisy measurements from a limited number of sensors.
#'
#' @param sites Integer. Number of independent soil sites.
#'
#' @param hours Numeric. Total simulation duration.
#'
#' @param depth Numeric. Maximum soil depth of the spatial domain.
#'
#' @param n_depth_points Integer. Number of spatial grid points used to
#' discretize the soil profile.
#'
#' @param mean_initial_c Numeric. Center location of the initial Gaussian
#' contaminant concentration profile.
#'
#' @param sd_initial_c Numeric. Width parameter of the initial Gaussian
#' concentration profile.
#'
#' @param mean_diffusivity Numeric. Mean diffusivity parameter before applying
#' site-level covariate effects.
#'
#' @param sd_sites_diffusivity Numeric. Standard deviation of unexplained
#' site-level diffusivity variation.
#'
#' @param beta_porosity Numeric. Effect of soil porosity on diffusivity.
#'
#' @param safety_factor Numeric between 0 and 1. Multiplier applied to the
#' maximum stable numerical time step determined from the diffusion CFL
#' condition.
#'
#' @param sensor_locations Numeric vector. Depth locations where observations
#' are collected.
#'
#' @param sensor_times Numeric vector. Times where observations are recorded.
#'
#' @param sigma_obs Numeric. Observation noise parameter.
#'
#' @return
#' A list containing:
#'
#' \describe{
#' \item{true_fields}{
#' Complete simulated PDE concentration fields for every site.
#' Contains the unobserved spatial process.
#' }
#'
#' \item{observations}{
#' Simulated sensor measurements collected only at selected depths and times.
#' }
#'
#' \item{parameters}{
#' Data frame containing site-level parameters including porosity and
#' diffusivity.
#' }
#' }
#'
#' @details
#' The simulation represents the causal structure:
#'
#' \deqn{
#' porosity_j \rightarrow D_j \rightarrow u_j(x,t)
#' \rightarrow y_{jkl}
#' }
#'
#' meaning soil structure influences contaminant diffusion through its effect
#' on the diffusion coefficient.
#'
#' The numerical solver automatically adjusts the integration time step based
#' on the diffusion stability condition:
#'
#' \deqn{
#' \frac{D\Delta t}{\Delta x^2}<0.5
#' }
#'
#' ensuring stable diffusion simulations.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_pde_contaminant_diffusion(
#'   sites = 15,
#'   hours = 72
#' )
#'
#' head(sim$observations)
#'
#' head(sim$parameters)
#' }
#'
#' @export
simulate_pde_contaminant_diffusion <- function(
  # Settings 
  sites = 15,
  hours = 48,
  depth = 3,
  n_depth_points = 30,

  ## Initial Field 
  mean_initial_c = 2,
  sd_initial_c   = 1,

  # Diffusivity
  mean_diffusivity     = 2,
  sd_sites_diffusivity = 1,
  beta_porosity        = 0.1,
  safety_factor        = 0.8,

  # Sensors
  sensor_locations = c(0.3, 1.05, 1.95, 2.7),
  sensor_times     = seq(0, hours, length.out = 12),
  sigma_obs        = 0.1
){

  # Simulate Covarite porosity
  porosity <- rnorm(sites, mean = 0, sd = 1)

  # Create a spatial grid 
  spatial_grid <- make_spatial_grid(length = depth, n_points = n_depth_points)

  # Create Initial Conditions 
  initial <- make_initial_gaussian(
    grid    = spatial_grid, 
    center  = mean_initial_c,
    width   = sd_initial_c
  )

  # Hierarchical diffusivity
  D <- log(mean_diffusivity) + beta_porosity * porosity + sd_sites_diffusivity * rnorm(sites,0, 1)
  D <- exp(D)

  # Set boundery conditions 
  boundary <- make_boundary_fixed(right = 0, left = 0) 
  
  # Simulate the PDE 
  results <- lapply(seq(sites), function(j){

    # Pick dt trying not to blow up 
    cfl <- check_cfl_diffusion(D = D[j], dx = spatial_grid$dx, dt = 1)
    dt_j <- safety_factor * cfl$max_stable_dt

    # Create a time grid with safe dt 
    time_grid <- make_time_grid(start = 0, end = hours, dt = dt_j)

    # PDE Equations 
    pde_eq <- make_diffusion_equation(D = D[j])

    # Solve the PDE 
    fields <- simulate_pde(
      grid = spatial_grid,
      time_grid = temporal_grid,
      initial_condition = initial,
      boundary_condition = boundary,
      equation = pde_eq,
      method = c("rk4"),
      record_every = max(1, round(1 / dt_j))
    )
    fields$site <- j

    # Observations 
    obs <- simulate_pde_observations(
      fields,
      sensor_locations = sensor_locations,
      sensor_times     = sensor_times,
      sigma_obs        = sigma_obs,
      likelihood       = "lognormal"
    )
    obs$site <- j

    # Return a list with true fields and observations 
    list(field = fields, obs = obs)
  }) 

  # Collect and return the results 
  list(
    true_fields  = do.call(rbind, lapply(results, `[[`, "field")),
    observations = do.call(rbind, lapply(results, `[[`, "obs")),
    parameters   = data.frame(site = seq_len(sites), porosity = porosity, D = D)
  )
    
}