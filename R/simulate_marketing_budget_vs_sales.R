

simulate_marketing_budget_vs_sales <- function(
  num_regions = 60,
  num_weeks_for_each_region = 30,
  prob_holiday_week         = 0.1,
  lambda_competion_density  = 3,
  budget_range              = c(1, 50),
  spline_ampliture          = 6,
  tau_budget                = 0.2,
  baseline_sales_per_week   = 20,
  sd_sales_region     = 3,
  baseline_comp_dens  = -1,
  sd_region_comp_desn = 0.3,
  beta_holiday        = -5,
  sales_sd            = 3
){
  # Simulate J number of regions each with N number of weekly sales
  ids <- make_nested_ids(levels = list(region = num_regions, num_weeks_for_each_region = num_weeks_for_each_region))
  
  # Simulate the Covariates Holiday indicator, Competion Density, Budget Spending 
  is_holiday        <- rbinom(n = nrow(ids), size = 1, prob = prob_holiday_week) 
  competion_density <- rpois(n  = nrow(ids), lambda = lambda_competion_density)
  budget            <- runif(n  = nrow(ids), min = budget_range[1], max = budget_range[2])

  # Build Non Linear Function P-Spine from budget 
  f_budget <- simulate_p_spline(x = budget, df = spline_ampliture, tau = tau_budget)

  # Build the Hierarcle paramters for the intercet and the competition density 
  alpha_j <- baseline_sales_per_week + sd_sales_region * rnorm(n = num_regions, mean = 0, sd = 1) 
  comp_j  <- baseline_comp_dens + sd_region_comp_desn  * rnorm(n = num_regions, mean = 0, sd = 1)

  # Build a Linear Predictor and Sample weekly sales from Normal Distribution 
  mu_i <- (
    alpha_j[ids$region_id] 
    + comp_j[ids$region_id] * competion_density
    + beta_holiday * is_holiday
    + f_budget$fitted
  )
  sales <- pmax(rnorm(n = nrow(ids), mean = mu_i, sd = sales_sd), 0)

  # Combine in one dataset 
  sim_data <- data.frame(
    region_id  = region,
    competion_density = competion_density,
    is_holiday        = is_holiday,
    budget            = budget,
    sales             = sales
  )
}