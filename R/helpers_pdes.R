# =============================================================================
# helpers_pde.R
#
# Composable "Lego block" helpers for simulating 1D PDEs, in the same spirit
# as helpers_state_space.R: grid / initial-condition / boundary-condition /
# equation are independent, swappable pieces, and one generic driver
# (simulate_pde) wires them together.
#
# Design changes vs. the previous version (and why):
#
#   1. make_spatial_grid() now actually tags its output with class
#      "pde_grid". Previously nothing set this class, so every
#      inherits(grid, "pde_grid") check downstream silently failed.
#
#   2. Boundary conditions are pure data (S3 objects: type + params), applied
#      through one generic, apply_boundary(). Previously boundary objects
#      were structs but simulate_pde() called boundary_condition(u, t) as if
#      they were functions -- that's a straight-up bug, not a style choice.
#
#   3. Equation constructors (make_diffusion_equation(), etc.) now return
#      ONLY a right-hand-side function du/dt = f(u, x, t, params). They used
#      to each hand-roll their own Euler step (state + lambda*(...)), which
#      duplicated the integrator inside every equation AND conflicted with
#      simulate_pde()'s own docstring, which already promised a generic
#      pde_rhs contract. Now there is exactly one place that does numerical
#      integration: the driver. Equations are just physics.
#
#   4. simulate_pde() itself had a scoping bug (used `grid` when the
#      argument was named `spatial_grid`) and no working relationship
#      between dt and time_grid. Rewritten as a single clean Euler/RK4
#      driver.
#
#   5. Added a CFL helper, because explicit diffusion schemes blow up
#      silently if dt is too large relative to dx^2/D, and that's a much
#      more useful failure mode to catch before simulating than after.
#
#   6. Second-order equations (the wave equation) genuinely don't fit the
#      first-order du/dt Lego block -- they need two state variables (u and
#      its velocity). Rather than force it into the same interface and hide
#      that mismatch, simulate_wave() is its own small driver that reuses
#      the grid/boundary blocks but is honest about needing extra state.
#
# =============================================================================

# ----------------------------------------------------------------------------
# 1. Grids
# ----------------------------------------------------------------------------

#' Create a Spatial Grid for PDE Simulations
#'
#' @description
#' Creates a regular 1D or 2D spatial grid used by PDE simulators. Returns an
#' object of class `pde_grid` storing coordinates plus spacing metadata.
#'
#' * For 1D: specify `length` and `n_points`.
#' * For 2D: specify `x_length`, `y_length`, `nx`, `ny`.
#'
#' @export
make_spatial_grid <- function(length = NULL, n_points = NULL,
                               x_length = NULL, y_length = NULL,
                               nx = NULL, ny = NULL) {

  if (!is.null(length) && !is.null(n_points)) {
    grid <- list(
      dimension = 1,
      x         = seq(0, length, length.out = n_points),
      dx        = length / (n_points - 1),
      length    = length,
      n_points  = n_points
    )
    class(grid) <- "pde_grid"
    return(grid)
  }

  if (!is.null(x_length) && !is.null(y_length) && !is.null(nx) && !is.null(ny)) {
    grid <- list(
      dimension = 2,
      x = seq(0, x_length, length.out = nx),
      y = seq(0, y_length, length.out = ny),
      dx = x_length / (nx - 1),
      dy = y_length / (ny - 1),
      x_length = x_length, y_length = y_length,
      nx = nx, ny = ny
    )
    class(grid) <- "pde_grid"
    return(grid)
  }

  stop("Supply either (length, n_points) for a 1D grid or ",
       "(x_length, y_length, nx, ny) for a 2D grid.")
}

#' Create a Temporal Simulation Grid
#'
#' @description
#' One-dimensional time axis for a PDE simulation. Provide exactly one of
#' `dt` (fixed step size) or `n_points` (fixed point count).
#'
#' @export
make_time_grid <- function(start = 0, end = NULL, dt = NULL, n_points = NULL) {
  if (!is.null(dt))       return(seq(start, end, by = dt))
  if (!is.null(n_points)) return(seq(start, end, length.out = n_points))
  stop("Provide either dt or n_points.")
}

#' Check the CFL Stability Condition for an Explicit Diffusion Step
#'
#' @description
#' Explicit (forward Euler) diffusion schemes are only numerically stable if
#' `dt <= dx^2 / (2*D)`. This never crashes on its own -- it just quietly
#' produces garbage (oscillation, blow-up) if violated -- so it's worth
#' checking before, not after, a long simulation run.
#'
#' @export
check_cfl_diffusion <- function(D, dx, dt) {
  max_stable_dt <- dx^2 / (2 * D)
  list(stable = dt <= max_stable_dt, dt = dt, max_stable_dt = max_stable_dt)
}

# ----------------------------------------------------------------------------
# 2. Initial conditions (1D)
# ----------------------------------------------------------------------------

#' Constant Initial Condition
#' @export
make_initial_constant <- function(grid, value = 0) {
  if (!inherits(grid, "pde_grid")) stop("grid must be created by make_spatial_grid().")
  if (grid$dimension == 1) return(rep(value, grid$n_points))
  matrix(value, nrow = grid$ny, ncol = grid$nx)
}

#' Point-Source Initial Condition
#' @export
make_initial_point_source <- function(grid, location, amplitude = 100, background = 0) {
  if (grid$dimension != 1) stop("Currently implemented only for 1D grids.")
  u0 <- rep(background, grid$n_points)
  idx <- which.min(abs(grid$x - location))
  u0[idx] <- amplitude
  u0
}

#' Gaussian Initial Condition
#' @export
make_initial_gaussian <- function(grid, center, width, amplitude = 1, background = 0) {
  if (grid$dimension != 1) stop("Currently implemented only for 1D grids.")
  background + amplitude * exp(-(grid$x - center)^2 / (2 * width^2))
}

#' Step Initial Condition
#' @export
make_initial_step <- function(grid, cutoff, left = 1, right = 0) {
  if (grid$dimension != 1) stop("Currently implemented only for 1D grids.")
  ifelse(grid$x < cutoff, left, right)
}

#' Sinusoidal Initial Condition
#' @export
make_initial_sine <- function(grid, amplitude = 1, wavelength, phase = 0) {
  if (grid$dimension != 1) stop("Currently implemented only for 1D grids.")
  amplitude * sin(2 * pi * grid$x / wavelength + phase)
}

#' Random (White-Noise) Initial Condition
#' @export
make_initial_random <- function(grid, mean = 0, sd = 1) {
  if (grid$dimension != 1) stop("Currently implemented only for 1D grids.")
  rnorm(grid$n_points, mean, sd)
}

#' Custom Initial Condition
#' @export
make_initial_custom <- function(grid, fun) {
  if (grid$dimension != 1) stop("Currently implemented only for 1D grids.")
  fun(grid$x)
}

# ----------------------------------------------------------------------------
# 3. Boundary conditions
# ----------------------------------------------------------------------------
# All boundary constructors return plain S3 data objects (class "pde_boundary").
# apply_boundary() is the single generic that knows how to enforce each type.
# This is the piece that was actually broken before: boundary objects were
# never callable, but simulate_pde() tried to call them.

#' Fixed (Dirichlet) Boundary Conditions
#' @export
make_boundary_fixed <- function(left = 0, right = 0) {
  structure(list(type = "fixed", left = left, right = right), class = "pde_boundary")
}

#' Reflecting (Zero-Flux / Neumann) Boundary Conditions
#' @export
make_boundary_reflecting <- function() {
  structure(list(type = "reflecting"), class = "pde_boundary")
}

#' Periodic Boundary Conditions
#' @export
make_boundary_periodic <- function() {
  structure(list(type = "periodic"), class = "pde_boundary")
}

#' Open (Outflow / Zero-Curvature) Boundary Conditions
#' @export
make_boundary_open <- function() {
  structure(list(type = "open"), class = "pde_boundary")
}

#' Custom Boundary Conditions
#'
#' @param fun function(u, t) returning the state after enforcing the
#'   boundary condition.
#' @export
make_boundary_custom <- function(fun) {
  structure(list(type = "custom", fun = fun), class = "pde_boundary")
}

#' Apply a Boundary Condition to a State Vector
#' @export
apply_boundary <- function(bc, u, t = NULL) {
  if (!inherits(bc, "pde_boundary")) {
    stop("bc must be created by one of the make_boundary_*() functions.")
  }
  n <- length(u)
  switch(bc$type,
    fixed      = { u[1] <- bc$left; u[n] <- bc$right; u },
    reflecting = { u[1] <- u[2];    u[n] <- u[n - 1];  u },
    periodic   = { u[1] <- u[n - 1]; u[n] <- u[2];     u },
    open       = { u[1] <- 2 * u[2] - u[3]; u[n] <- 2 * u[n - 1] - u[n - 2]; u },
    custom     = bc$fun(u, t),
    stop("Unknown boundary type: ", bc$type)
  )
}

# ----------------------------------------------------------------------------
# 4. Equations: each returns a pde_equation whose $rhs computes du/dt only
# ----------------------------------------------------------------------------

#' @keywords internal
new_pde_equation <- function(name, rhs, parameters = list()) {
  structure(list(name = name, rhs = rhs, parameters = parameters), class = "pde_equation")
}

#' @keywords internal
validate_pde_equation <- function(equation) {
  if (!inherits(equation, "pde_equation")) stop("Object must inherit from class 'pde_equation'.")
  if (!is.function(equation$rhs)) stop("PDE equation requires an rhs function.")
  invisible(TRUE)
}

#' Second Spatial Derivative (Laplacian), Central Differences, Vectorised
#'
#' @description
#' Edge entries are filled with their nearest interior neighbour purely as a
#' placeholder -- they get overwritten by apply_boundary() immediately after,
#' so their exact value here never matters for a correctly-configured
#' simulation.
#' @export
laplacian_1d <- function(u, dx) {
  n <- length(u)
  lap <- numeric(n)
  lap[2:(n - 1)] <- (u[1:(n - 2)] - 2 * u[2:(n - 1)] + u[3:n]) / dx^2
  lap[1] <- lap[2]
  lap[n] <- lap[n - 1]
  lap
}

#' First Spatial Derivative, Upwind Scheme
#'
#' @description
#' Direction is chosen from the sign of `v` (single scalar velocity, as in
#' the site-level advection models below): backward difference for v >= 0
#' (transport moving in the +x direction), forward difference for v < 0.
#' @export
upwind_grad <- function(u, dx, v) {
  n <- length(u)
  grad <- numeric(n)
  if (v >= 0) {
    grad[2:n] <- (u[2:n] - u[1:(n - 1)]) / dx
    grad[1] <- grad[2]
  } else {
    grad[1:(n - 1)] <- (u[2:n] - u[1:(n - 1)]) / dx
    grad[n] <- grad[n - 1]
  }
  grad
}

#' Diffusion (Heat) Equation: du/dt = D * d2u/dx2
#' @export
make_diffusion_equation <- function(D = 1) {
  rhs <- function(u, grid, t, parameters) {
    parameters$D * laplacian_1d(u, grid$dx)
  }
  new_pde_equation("diffusion", rhs, list(D = D))
}

#' Pure Advection (Transport) Equation: du/dt = -v * du/dx
#' @export
make_advection_equation <- function(v = 1) {
  rhs <- function(u, grid, t, parameters) {
    -parameters$v * upwind_grad(u, grid$dx, parameters$v)
  }
  new_pde_equation("advection", rhs, list(v = v))
}

#' Advection-Diffusion Equation: du/dt = -v*du/dx + D*d2u/dx2
#' @export
make_advection_diffusion_equation <- function(v = 1, D = 1) {
  rhs <- function(u, grid, t, parameters) {
    v <- parameters$v; D <- parameters$D
    -v * upwind_grad(u, grid$dx, v) + D * laplacian_1d(u, grid$dx)
  }
  new_pde_equation("advection_diffusion", rhs, list(v = v, D = D))
}

#' Reaction-Diffusion (Fisher-KPP / Logistic Growth) Equation:
#' du/dt = D*d2u/dx2 + r*u*(1 - u/K)
#' @export
make_reaction_diffusion_equation <- function(D = 1, r = 1, K = 1) {
  rhs <- function(u, grid, t, parameters) {
    D <- parameters$D; r <- parameters$r; K <- parameters$K
    D * laplacian_1d(u, grid$dx) + r * u * (1 - u / K)
  }
  new_pde_equation("reaction_diffusion", rhs, list(D = D, r = r, K = K))
}

#' Custom Equation
#'
#' @param rhs function(u, grid, t, parameters) returning du/dt.
#' @export
make_custom_equation <- function(rhs, parameters = list(), name = "custom") {
  new_pde_equation(name, rhs, parameters)
}

# ----------------------------------------------------------------------------
# 5. The generic 1st-order driver
# ----------------------------------------------------------------------------

#' Simulate a 1D First-Order PDE
#'
#' @description
#' Generic driver: takes a grid, a time grid, an initial condition, a
#' boundary condition, and an equation (an rhs function bundle) and produces
#' a long data frame of (time, x, value). This is the one place time
#' integration happens; equations only ever describe du/dt.
#'
#' @param grid A `pde_grid` from make_spatial_grid().
#' @param time_grid A numeric vector from make_time_grid().
#' @param initial_condition Numeric vector, or function(x) -> numeric vector.
#' @param boundary_condition A `pde_boundary` object.
#' @param equation A `pde_equation` object.
#' @param method "euler" (fast, needs CFL-respecting dt) or "rk4" (more
#'   forgiving of larger dt, ~4x the cost per step).
#' @param record_every Save output every this many steps.
#'
#' @export
simulate_pde <- function(grid, time_grid, initial_condition, boundary_condition,
                          equation, method = c("euler", "rk4"), record_every = 1) {

  method <- match.arg(method)
  validate_pde_equation(equation)

  u <- if (is.function(initial_condition)) initial_condition(grid$x) else initial_condition
  if (length(u) != length(grid$x)) {
    stop("Initial condition length does not match spatial grid.")
  }

  dt <- diff(time_grid)[1]
  if (length(time_grid) > 2 && any(abs(diff(diff(time_grid))) > 1e-8 * dt)) {
    warning("time_grid is not evenly spaced; using dt from the first interval throughout.")
  }

  f <- function(u, t) equation$rhs(u, grid, t, equation$parameters)

  n_steps <- length(time_grid)
  output <- vector("list", n_steps)
  counter <- 0L

  for (step in seq_len(n_steps)) {
    t <- time_grid[step]

    if ((step - 1) %% record_every == 0) {
      counter <- counter + 1L
      output[[counter]] <- data.frame(time = t, x = grid$x, value = u)
    }

    if (step == n_steps) break

    u_new <- switch(method,
      euler = u + dt * f(u, t),
      rk4   = {
        k1 <- f(u, t)
        k2 <- f(u + dt / 2 * k1, t + dt / 2)
        k3 <- f(u + dt / 2 * k2, t + dt / 2)
        k4 <- f(u + dt * k3, t + dt)
        u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      }
    )

    u <- apply_boundary(boundary_condition, u_new, t)
  }

  do.call(rbind, output[seq_len(counter)])
}

# ----------------------------------------------------------------------------
# 6. Wave equation: its own small driver (genuinely second-order)
# ----------------------------------------------------------------------------

#' Simulate the 1D Wave Equation: d2u/dt2 = c^2 * d2u/dx2
#'
#' @description
#' Kept separate from simulate_pde() on purpose: the wave equation needs two
#' state variables (displacement u and velocity du/dt), which doesn't fit the
#' single-vector du/dt contract the other equations share. Forcing it into
#' the same interface would mean hiding a real structural difference instead
#' of documenting it.
#'
#' @param initial_v Initial velocity field (defaults to zero, i.e. "released
#'   from rest").
#' @export
simulate_wave <- function(grid, time_grid, initial_u, initial_v = NULL,
                           boundary_condition, c = 1, record_every = 1) {

  if (is.null(initial_v)) initial_v <- rep(0, grid$n_points)

  u <- if (is.function(initial_u)) initial_u(grid$x) else initial_u
  v <- initial_v
  dt <- diff(time_grid)[1]

  n_steps <- length(time_grid)
  output <- vector("list", n_steps)
  counter <- 0L

  for (step in seq_len(n_steps)) {
    t <- time_grid[step]

    if ((step - 1) %% record_every == 0) {
      counter <- counter + 1L
      output[[counter]] <- data.frame(time = t, x = grid$x, value = u)
    }

    if (step == n_steps) break

    accel <- c^2 * laplacian_1d(u, grid$dx)
    u_new <- u + dt * v
    v      <- v + dt * accel
    u      <- apply_boundary(boundary_condition, u_new, t)
  }

  do.call(rbind, output[seq_len(counter)])
}

# ----------------------------------------------------------------------------
# 7. Observation model: PDE field -> sparse, noisy sensor readings
# ----------------------------------------------------------------------------

#' Generate Observations from a PDE Solution
#'
#' @description
#' Converts a simulated PDE field into noisy, sparsely-sampled observations,
#' i.e. the layer between the PDE simulator and a statistical inference
#' model. Supports both a generic random/regular-grid subsample (as before)
#' and, importantly for sensor-network-style data, snapping to a fixed set
#' of `sensor_locations` (and optionally `sensor_times`) -- e.g. "the 3-4
#' depths where a probe actually sits", which the old grid/random sampling
#' modes couldn't express directly.
#'
#' @param likelihood "normal" (additive Gaussian noise) or "lognormal"
#'   (multiplicative, always-positive noise -- appropriate for
#'   concentrations/counts that can't go negative).
#' @param sensor_locations Optional vector of spatial coordinates; each is
#'   snapped to the nearest grid location actually present in
#'   `pde_solution`.
#' @param sensor_times Optional vector of times, snapped the same way. If
#'   omitted, all simulated times are kept for the chosen locations.
#'
#' @export
simulate_pde_observations <- function(pde_solution,
                                       num_locations = 5, num_times = 20,
                                       sigma_obs = 0.1,
                                       sampling_type = c("random", "grid"),
                                       likelihood = c("normal", "lognormal"),
                                       sensor_locations = NULL,
                                       sensor_times = NULL) {

  sampling_type <- match.arg(sampling_type)
  likelihood <- match.arg(likelihood)

  if (!all(c("time", "x", "value") %in% names(pde_solution))) {
    stop("pde_solution must contain time, x, and value columns.")
  }

  if (!is.null(sensor_locations)) {
    ux <- unique(pde_solution$x)
    snapped_x <- vapply(sensor_locations, function(loc) ux[which.min(abs(ux - loc))], numeric(1))
    sampled <- pde_solution[pde_solution$x %in% snapped_x, ]

    if (!is.null(sensor_times)) {
      ut <- unique(pde_solution$time)
      snapped_t <- vapply(sensor_times, function(tt) ut[which.min(abs(ut - tt))], numeric(1))
      sampled <- sampled[sampled$time %in% snapped_t, ]
    }

  } else if (sampling_type == "random") {
    idx <- sample(seq_len(nrow(pde_solution)),
                   size = min(num_locations * num_times, nrow(pde_solution)))
    sampled <- pde_solution[idx, ]

  } else { # "grid"
    ux <- unique(pde_solution$x)
    ut <- unique(pde_solution$time)
    x_sel <- ux[round(seq(1, length(ux), length.out = num_locations))]
    t_sel <- ut[round(seq(1, length(ut), length.out = num_times))]
    sampled <- pde_solution[pde_solution$x %in% x_sel & pde_solution$time %in% t_sel, ]
  }

  sampled$true_value <- sampled$value

  sampled$observation <- if (likelihood == "normal") {
    rnorm(nrow(sampled), mean = sampled$value, sd = sigma_obs)
  } else {
    rlnorm(nrow(sampled), meanlog = log(pmax(sampled$value, 1e-8)), sdlog = sigma_obs)
  }

  sampled[, c("time", "x", "true_value", "observation")]
}