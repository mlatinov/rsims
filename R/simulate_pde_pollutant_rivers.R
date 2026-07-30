simulate_pde_pollutant_rivers <- function(
  # Settings 
  rivers = 15,
  river_length  = 26,
  n_depth_points = 20,
  hours  = 50,
  mean_initial = 0.3,
  sd_initial   = 0.1,
  safety_factor = 0.8,

  # Covariates 
  mean_flow_rate = 1,
  mean_channel_roughness = 2,

  # Covariates Paramters effects over PDE 
  mean_rate_v = -0.2,
  beta_flow   = 0.3,
  sigma_rivers_v = 0.1,
  mean_diffusion = -2.3,
  beta_roughness = 0.1,
  sigma_rivers_diff = 0.1,

  # Stations 
  station_location = c(1, 5, 10, 15, 20),
  station_times    = seq(0, hours, length.out=12),

  # Observations
  sigma_obs = 0.1

){

  # Simulate Covariates 
  flow_rate         <- rexp(rivers, rate = mean_flow_rate)
  channel_roughness <- rlnorm(rivers, meanlog = mean_channel_roughness, sd = 0.5)

  # Create a spatical location grid 
  location_grid <- make_spatial_grid(length = river_length, n_points = n_depth_points)

  # Create the initial condutions 
  initial <- make_initial_gaussian(grid = location_grid, center = mean_initial, width = sd_initial)

  # Hierarchical Structure of the PDE paramters 
  D <- exp(mean_diffusion + beta_roughness * log(channel_roughness) + sigma_rivers_diff * rnorm(rivers, 0, 1))
  v <- exp(mean_rate_v + beta_flow * log(flow_rate) + sigma_rivers_v * rnorm(rivers, 0, 1)) 

  # Set River Boundries 
  boundary <- make_boundary_fixed(left = 0, right = 0)

  # Simulate PDE 
  results <- lapply(seq(rivers), function(j){

    # Pick safe dt 
    cfl <- check_cfl_diffusion(D = D[j],dx = location_grid$dx, dt = 1)
    dt_j <- safety_factor * cfl$max_stable_dt

    # Create a temporal grid with safe dt
    temporal_grid <- make_time_grid(start = 0, end = hours, dt = dt_j)

    # Equation : Advection–diffusion
    pde_eq <- make_advection_diffusion_equation(v = v[j], D = D[j])

    # Solve for River J
    field <- simulate_pde(
      grid = location_grid,
      time_grid = temporal_grid,
      initial_condition = initial,
      boundary_condition = boundary,
      equation = pde_eq,
      method = "rk4",
      record_every = max(1, round(1 / dt_j))
    )
    field$river <- j

    # Generate Observations 
    obs <- simulate_pde_observations(
      pde_solution = field,
      sensor_locations =  station_location,
      sensor_times     =  station_times,
      sigma_obs        =  sigma_obs,
      likelihood       = "lognormal" 
    )
    obs$river <- j

    # Return a list with true fields and observations 
    list(field = field, obs = obs)
  })

  # Return a list with the collected true and observed fields 
  list(
    true_fields  = do.call(rbind, lapply(results, `[[`, "field")),
    observations = do.call(rbind, lapply(results, `[[`, "obs")),
    parameters   = data.frame(
      site = seq_len(rivers), 
      flow_rate = flow_rate,
      channel_roughness = channel_roughness, 
      D = D,
      v = v
    )
  )
}
