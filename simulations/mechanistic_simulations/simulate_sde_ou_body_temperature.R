simulate_sde_ou_body_temperature <- function(
    num_subjects       = 40,
    num_hours          = 48,
    mean_ambient_t = 20,                sd_ambient_t   = 4,
    mean_activity_l  = 5,               sd_activity_l    = 1.5,
    mean_initial_point = 36.8,          sd_initial_point = 0.5,
    mean_outflow_rate = log(0.15),   
    sd_subjects_outflow_rate = 0.08,
    beta_activity = -0.03,
    mean_heat_inflow_rate = log(0.15) + log(37), 
    sd_subjects_heat_inflow_rate = 0.08,
    beta_ambient_t = 0.015,
    body_temp_sd = 0.15,
    dt       = 1,
    ou_sigma = 0.1 
){
  
  # Simulate Subjects measurment of body temperature every hour to max hours
  n   <- num_subjects * num_hours 
  ids <- make_nested_ids(levels = list(subject = num_subjects, hour = num_hours)) 

  # Covariates 
  ambient_temp   <- rnorm(num_subjects, mean = mean_ambient_t, sd = sd_ambient_t)
  activity_level <- rnorm(num_subjects, mean = mean_activity_l, sd = sd_activity_l)

  # Center the covariates 
  ambient_c  <- ambient_temp   - mean_ambient_t
  activity_c <- activity_level - mean_activity_l

  # Paramters 
  a  <- exp(mean_heat_inflow_rate + beta_ambient_t * ambient_c + sd_subjects_heat_inflow_rate * rnorm(num_subjects, 0, 1))
  b  <- exp(mean_outflow_rate     + beta_activity  * activity_c + sd_subjects_outflow_rate * rnorm(num_subjects, 0 , 1))
  xo <- rnorm(num_subjects, mean = mean_initial_point, sd = sd_initial_point)

  # Simulate OU Process
  C <- matrix(NA_real_, nrow = num_hours, ncol = num_subjects)
  
  # Initilize the process 
  C[1, ] <- xo

  # Continue the OU Process
  for (i in 2:num_hours) {
    mean_t  <- a / b + (C[i - 1, ] - a / b) * exp(-b * dt)
    sigma_t <- sqrt((ou_sigma^2 / (2 * b)) * (1 - exp(-2 * b * dt)))
    C[i, ]  <- rnorm(num_subjects, mean = mean_t, sd = sigma_t)
  } 

  # Flatten the matrix and Sample the observations 
  latent_body_temp <- as.vector(C)
  body_temp        <- rnorm(n, mean = latent_body_temp, sd = body_temp_sd)

  # Return a dataframe with the simulated data 
  data.frame(
    subject  = ids$subject_id,
    hour     = ids$hour_id,
    ambient_t = ambient_temp[ids$subject_id],
    activity_level = activity_level[ids$subject_id],
    body_temp      = body_temp
  )
  
}
  