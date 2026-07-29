# ---------------------------------------------------------------------------
# 0. Core primitive: everything else is built on this
# ---------------------------------------------------------------------------

#' Simulate a Gaussian random walk
#'
#' The single most reused primitive in the whole file. A local level, a
#' trend's slope, a time-varying regression coefficient, and a spline
#' random walk are all just calls to this function.
#'
#' x_t = x_{t-1} + drift + eps_t,   eps_t ~ N(0, sigma^2)
#'
#' @param T integer, number of time points to return (t = 1..T)
#' @param sigma innovation standard deviation (scalar or length-T vector for
#'   time-varying volatility)
#' @param initial starting value x_0
#' @param drift constant drift added at every step (default 0)
#'
#' @return numeric vector of length T
#' @export
simulate_random_walk <- function(T, sigma = 1, initial = 0, drift = 0) {
  stopifnot(T >= 1)
  eps <- stats::rnorm(T, mean = 0, sd = sigma)
  initial + drift * seq_len(T) + cumsum(eps)
}


# ---------------------------------------------------------------------------
# 1. simulate_local_level()
# ---------------------------------------------------------------------------

#' Simulate a local level (random walk) state
#'
#' mu_t = mu_{t-1} + eps_t
#'
#' @param T integer horizon
#' @param initial starting level mu_0
#' @param sigma innovation sd
#'
#' @return list(state = mu, innovations = epsilon)
#' @examples
#' level <- simulate_local_level(T = 365, initial = log(50), sigma = 0.02)
#' @export
simulate_local_level <- function(T, initial = 0, sigma = 1) {
  epsilon <- stats::rnorm(T, mean = 0, sd = sigma)
  mu <- initial + cumsum(epsilon)
  list(state = mu, innovations = epsilon)
}


# ---------------------------------------------------------------------------
# 2. simulate_local_trend()
# ---------------------------------------------------------------------------

#' Simulate a local linear trend state
#'
#' mu_t = mu_{t-1} + b_{t-1} + eps_t
#' b_t  = b_{t-1} + eta_t
#'
#' Both recursions are vectorized: b is a random walk on the slope, and mu
#' is the cumulative sum of the *lagged* slope plus its own innovations.
#'
#' @param T integer horizon
#' @param level0 starting level mu_0
#' @param slope0 starting slope b_0
#' @param sigma_level innovation sd for the level equation
#' @param sigma_slope innovation sd for the slope equation
#'
#' @return list(level = mu, slope = b)
#' @examples
#' trend <- simulate_local_trend(T = 365, level0 = 100, slope0 = 0.3,
#'                                sigma_level = 1, sigma_slope = 0.05)
#' @export
simulate_local_trend <- function(T, level0 = 0, slope0 = 0,
                                  sigma_level = 1, sigma_slope = 1) {
  eta <- stats::rnorm(T, mean = 0, sd = sigma_slope)
  b <- slope0 + cumsum(eta)                 # b_1 .. b_T
  b_lag <- c(slope0, b[-T])                 # b_0 .. b_{T-1}

  eps <- stats::rnorm(T, mean = 0, sd = sigma_level)
  mu <- level0 + cumsum(b_lag + eps)

  list(level = mu, slope = b)
}


# ---------------------------------------------------------------------------
# 3. simulate_seasonality()
# ---------------------------------------------------------------------------

#' Recycle a seasonal effects pattern out to an arbitrary horizon
#'
#' Wraps the `rep(effects, length.out = T)` pattern you'd otherwise retype
#' in every script, plus optional centering so the seasonal component sums
#' to (approximately) zero and doesn't fight with the level/intercept.
#'
#' @param period length of one seasonal cycle (e.g. 7 for day-of-week)
#' @param effects numeric vector of length `period` with the effect for
#'   each phase of the cycle
#' @param T horizon to recycle out to. Defaults to `period` (one full
#'   cycle) if not supplied.
#' @param center if TRUE, subtract mean(effects) so the pattern is centered
#'   before recycling
#'
#' @return numeric vector of length T (named `season` when used in `$`
#'   contexts downstream)
#' @examples
#' weekly <- simulate_seasonality(period = 7,
#'                                 effects = c(-.1, -.2, -.1, 0, .15, .25, .10),
#'                                 T = 365, center = TRUE)
#' @export
simulate_seasonality <- function(period, effects, T = NULL, center = TRUE) {
  stopifnot(length(effects) == period)
  if (center) effects <- effects - mean(effects)
  if (is.null(T)) T <- period
  rep(effects, length.out = T)
}

# ---------------------------------------------------------------------------
# 3b. simulate_regime() -- the block-constant sibling of simulate_seasonality()
# ---------------------------------------------------------------------------
 
#' Simulate a piecewise-constant regime / structural-break series
#'
#' The mirror image of simulate_seasonality(): instead of tiling a short
#' repeating pattern many times across T (fine-grained, `rep(effects,
#' length.out = T)`), this HOLDS each value constant for `block_length`
#' periods before switching to the next one (coarse-grained, `rep(effects,
#' each = block_length, length.out = T)`).
#'
#' Use simulate_seasonality() when every unit *within* one cycle differs
#' (e.g. Monday vs Tuesday vs Wednesday, repeated every week).
#' Use simulate_regime() when a single value holds for a whole block, and
#' blocks change less often than every observation (e.g. a reservoir
#' sitting at one operating level for 50 days, then switching).
#'
#' Number of distinct effects needed = T / block_length (one per block).
#' If that isn't a whole number, the final block is truncated and a
#' warning is issued so it never fails silently.
#'
#' @param T integer horizon
#' @param block_length length of each regime/block in time units
#' @param effects numeric vector, one value per block. Length should equal
#'   `T / block_length`; if there are fewer effects than blocks, they
#'   recycle (e.g. ABAB...); if there are more, the extras are unused.
#' @param center if TRUE, subtract mean(effects) so the pattern is centered
#'   before expanding
#'
#' @return numeric vector of length T
#' @examples
#' # 4 regimes of 50 days each across a 200-day horizon
#' regime <- simulate_regime(T = 200, block_length = 50,
#'                            effects = c(0, -4, 5, 3))
#' @export
simulate_regime <- function(T, block_length, effects, center = FALSE) {
  n_blocks <- T / block_length
  if (n_blocks != round(n_blocks)) {
    warning(sprintf(
      "T (%d) is not a multiple of block_length (%d); final block will be truncated.",
      T, block_length
    ))
  }
  if (center) effects <- effects - mean(effects)
  rep(effects, each = block_length, length.out = T)
}
# ---------------------------------------------------------------------------
# 4. simulate_intervention()
# ---------------------------------------------------------------------------

#' Simulate an intervention / event indicator series
#'
#' Replaces the one-off `as.numeric(time >= t0)` (and friends) that shows
#' up in nearly every causal-impact-style simulation.
#'
#' @param T integer horizon
#' @param time the index at which the intervention starts (1-based)
#' @param type one of "step", "pulse", "ramp", "decay"
#' @param duration for "pulse": how many periods it stays on (default 1).
#'   unused for "ramp" (which always increments 1, 2, 3, ... to the end
#'   of the series) and "decay".
#' @param decay_rate for "decay": multiplicative decay per period (0 < r < 1)
#' @param magnitude scalar multiplier applied to the whole series (default 1)
#'
#' @return numeric vector of length T
#' @examples
#' step  <- simulate_intervention(T = 14, time = 8, type = "step")
#' pulse <- simulate_intervention(T = 14, time = 6, type = "pulse", duration = 1)
#' ramp  <- simulate_intervention(T = 14, time = 6, type = "ramp")
#' decay <- simulate_intervention(T = 14, time = 6, type = "decay", decay_rate = 0.8)
#' @export
simulate_intervention <- function(T, time, type = c("step", "pulse", "ramp", "decay"),
                                   duration = NULL, decay_rate = 0.8, magnitude = 1) {
  type <- match.arg(type)
  stopifnot(time >= 1, time <= T)
  t_idx <- seq_len(T)
  on <- t_idx >= time

  x <- switch(type,
    step = as.numeric(on),

    pulse = {
      dur <- if (is.null(duration)) 1 else duration
      as.numeric(on & t_idx < time + dur)
    },

    ramp = {
      # integer ramp: 0 0 0 0 0 1 2 3 4 5 ... (matches the spec's example)
      out <- numeric(T)
      out[on] <- seq_len(sum(on))
      out
    },

    decay = {
      out <- numeric(T)
      idx <- t_idx[on]
      out[on] <- decay_rate ^ (idx - time)
      out
    }
  )

  magnitude * x
}


# ---------------------------------------------------------------------------
# 5. simulate_state_regression()
# ---------------------------------------------------------------------------

#' Add a regression component onto (optionally) an existing state
#'
#' eta = state + X %*% betas
#'
#' @param state optional numeric vector to add the regression term to
#'   (e.g. a level or trend). If NULL, only X %*% betas is returned.
#' @param X a data.frame or matrix of predictors, one column per beta
#' @param betas numeric vector of coefficients, one per column of X
#'
#' @return numeric vector `eta`, same length as nrow(X)
#' @examples
#' eta <- simulate_state_regression(
#'   state = level$state,
#'   X = data.frame(temp = rnorm(365), holiday = rbinom(365, 1, .05)),
#'   betas = c(-0.01, 0.25)
#' )
#' @export
simulate_state_regression <- function(state = NULL, X, betas) {
  X <- as.matrix(X)
  stopifnot(ncol(X) == length(betas))
  reg <- as.numeric(X %*% betas)
  if (is.null(state)) return(reg)
  stopifnot(length(state) == length(reg))
  state + reg
}


# ---------------------------------------------------------------------------
# 6. simulate_bsts() -- the assembly point
# ---------------------------------------------------------------------------

#' Assemble state + seasonality + regression + intervention into observed data
#'
#' Sums whatever components you give it on the linear-predictor scale, then
#' draws from the requested observation family. Any component you don't
#' need can simply be omitted (defaults to 0).
#'
#' @param T integer horizon (required so the function knows how long to
#'   simulate if all components are scalars/NULL)
#' @param state numeric vector, e.g. `$state` from simulate_local_level() or
#'   `$level` from simulate_local_trend()
#' @param seasonal numeric vector from simulate_seasonality()
#' @param regression numeric vector from simulate_state_regression()
#' @param intervention numeric vector from simulate_intervention()
#' @param family one of "gaussian", "poisson", "binomial"
#' @param sigma_obs observation noise sd, only used when family = "gaussian"
#'
#' @return list(y = observed series, eta = linear predictor,
#'              components = the list of inputs actually used)
#' @examples
#' lvl <- simulate_local_level(T = 365, initial = log(50), sigma = 0.02)
#' wk  <- simulate_seasonality(period = 7,
#'                              effects = c(-.1,-.2,-.1,0,.15,.25,.10),
#'                              T = 365)
#' sim <- simulate_bsts(T = 365, state = lvl$state, seasonal = wk,
#'                       family = "poisson")
#' @export
simulate_bsts <- function(T,
                           state = NULL,
                           seasonal = NULL,
                           regression = NULL,
                           intervention = NULL,
                           family = c("gaussian", "poisson", "binomial"),
                           sigma_obs = 1) {
  family <- match.arg(family)

  zero_if_null <- function(x) if (is.null(x)) 0 else x

  eta <- zero_if_null(state) + zero_if_null(seasonal) +
    zero_if_null(regression) + zero_if_null(intervention)

  # broadcast scalars up to length T if every component was NULL/scalar
  if (length(eta) == 1L) eta <- rep(eta, T)

  y <- switch(family,
    gaussian = stats::rnorm(length(eta), mean = eta, sd = sigma_obs),
    poisson  = stats::rpois(length(eta), lambda = exp(eta)),
    binomial = stats::rbinom(length(eta), size = 1, prob = stats::plogis(eta))
  )

  list(
    y = y,
    eta = eta,
    components = list(
      state = state, seasonal = seasonal,
      regression = regression, intervention = intervention
    )
  )
}


# ---------------------------------------------------------------------------
# 7. simulate_hierarchical_states()
# ---------------------------------------------------------------------------

#' Simulate independent local-level states for many groups (e.g. stores)
#'
#' Each group j gets its own random-walk state, with group-level starting
#' points drawn from N(initial_mean, sigma_between^2) and within-group
#' innovation variance sigma_state^2.
#'
#' state j, t = mu_j0 + cumsum(eps_{j,t}),  eps_{j,t} ~ N(0, sigma_state^2)
#' mu_j0 ~ N(initial_mean, sigma_between^2)
#'
#' @param groups integer number of groups (e.g. 30 stores)
#' @param T integer horizon
#' @param initial_mean grand mean of the group-level starting points
#' @param sigma_between sd of group-level starting points around
#'   initial_mean
#' @param sigma_state innovation sd of each group's local level
#'
#' @return a `groups x T` numeric matrix, `state[group, time]`
#' @examples
#' states <- simulate_hierarchical_states(groups = 30, T = 365,
#'                                         initial_mean = log(50),
#'                                         sigma_between = .2,
#'                                         sigma_state = .01)
#' @export
simulate_hierarchical_states <- function(groups, T, initial_mean = 0,
                                          sigma_between = 1, sigma_state = 1) {
  group_starts <- stats::rnorm(groups, mean = initial_mean, sd = sigma_between)
  eps <- matrix(stats::rnorm(groups * T, mean = 0, sd = sigma_state),
                nrow = groups, ncol = T)
  state <- group_starts + t(apply(eps, 1, cumsum))
  rownames(state) <- paste0("group_", seq_len(groups))
  state
}


# ---------------------------------------------------------------------------
# 8. simulate_multistate() -- correlated state evolution
# ---------------------------------------------------------------------------

#' Simulate correlated multivariate random-walk states
#'
#' Generalizes simulate_local_trend()-style coupling to k states whose
#' innovations are jointly Gaussian, e.g. level & slope moving together,
#' or several correlated regime states.
#'
#' x_t = x_{t-1} + MVN(0, Sigma)
#'
#' @param T integer horizon
#' @param initial numeric vector of length k, starting values
#' @param Sigma k x k innovation covariance matrix
#' @param drift optional numeric vector of length k, constant drift added
#'   at every step (default: no drift)
#'
#' @return `T x k` numeric matrix, one column per state
#' @examples
#' # variances 1 and 0.25 (sds 1 and 0.5); covariance 0.3 => corr = 0.6
#' Sigma <- matrix(c(1, 0.3, 0.3, 0.25), nrow = 2)
#' xs <- simulate_multistate(T = 200, initial = c(100, 0.3), Sigma = Sigma)
#' @export
simulate_multistate <- function(T, initial, Sigma, drift = NULL) {
  k <- length(initial)
  stopifnot(all(dim(Sigma) == c(k, k)))
  if (is.null(drift)) drift <- rep(0, k)

  # Cholesky factor so we don't need a MASS::mvrnorm dependency
  L <- chol(Sigma)                       # upper-triangular: Sigma = t(L) %*% L
  Z <- matrix(stats::rnorm(T * k), nrow = T, ncol = k)
  innovations <- Z %*% L                 # T x k, correlated innovations

  drift_mat <- matrix(drift, nrow = T, ncol = k, byrow = TRUE)
  increments <- drift_mat + innovations
  state <- sweep(apply(increments, 2, cumsum), 2, initial, `+`)

  colnames(state) <- paste0("state_", seq_len(k))
  state
}