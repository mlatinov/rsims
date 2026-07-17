## ============================================================================
##  Corrected spline-effect simulation helpers
##
##  KEY DESIGN DECISION (applies to every function):
##  The returned `weights` and `fitted` must satisfy  fitted == B %*% weights
##  EXACTLY. The old code rescaled `fitted` (f/sd(f)*amplitude) but returned the
##  UN-scaled weights, so the returned "truth" did not reproduce the output --
##  fatal for a simulate -> fit -> recover workflow.
##
##  FIX: when a target amplitude is requested, we scale the WEIGHTS, then compute
##  fitted from the scaled weights. Invariant f == B %*% w always holds.
##
##  A safe-sd() helper guards the divide-by-zero when a curve is near-constant.
## ============================================================================

## internal: standard deviation with a floor, never returns 0
.safe_sd <- function(v) {
  s <- stats::sd(v)
  if (!is.finite(s) || s < 1e-8) 1 else s
}

## internal: build a first-order random-walk vector, non-centered
##   w[1] = tau*z[1];  w[k] = w[k-1] + tau*z[k]
.rw1 <- function(K, tau) {
  z <- stats::rnorm(K)
  w <- numeric(K)
  w[1] <- tau * z[1]
  if (K > 1) for (k in 2:K) w[k] <- w[k - 1] + tau * z[k]
  w
}


#' Simulate a Smooth Effect from a Basis Matrix
#'
#' @description
#' Generates a smooth nonlinear effect by combining a basis matrix with spline
#' coefficients. By default the coefficients follow a first-order random walk,
#' producing a random smooth function; alternatively supply your own `weights`
#' to reuse a fixed nonlinear function across simulation studies.
#'
#' Random-walk coefficients satisfy
#' \deqn{w_1 = \tau z_1, \qquad w_k = w_{k-1} + \tau z_k, \quad z_k \sim N(0,1).}
#' The effect is \eqn{f(x) = Bw}.
#'
#' If `amplitude` is not `NULL`, the **weights are scaled** so that
#' \eqn{sd(f) = }`amplitude`. Because the weights (not the output) are scaled,
#' the returned `weights` and `fitted` always satisfy `fitted == B %*% weights`,
#' so the returned coefficients are the true generating coefficients -- suitable
#' for recovery checks.
#'
#' @param B Numeric basis matrix (rows = observations, columns = basis funs).
#' @param tau Positive random-walk step size (smoothness). Used only when
#'   `weights = NULL`.
#' @param amplitude Target standard deviation of the effect, achieved by scaling
#'   the weights. Set to `NULL` to leave the raw random-walk scale untouched.
#' @param weights Optional coefficient vector of length `ncol(B)`. If supplied,
#'   used directly (and still scaled to `amplitude` unless `amplitude = NULL`).
#'
#' @return List with `basis`, `weights` (the TRUE generating coefficients), and
#'   `fitted` (equal to `B %*% weights`).
#'
#' @examples
#' x <- sort(runif(200)); B <- splines::bs(x, df = 8)
#' s <- simulate_spline_effect(B)
#' stopifnot(all.equal(as.numeric(B %*% s$weights), s$fitted))  # invariant holds
#' @export
simulate_spline_effect <- function(
    B,
    tau       = 0.10,
    amplitude = 5,
    weights   = NULL
){
  K <- ncol(B)

  if (is.null(weights)) {
    w <- .rw1(K, tau)
  } else {
    if (length(weights) != K)
      stop("'weights' must have length equal to ncol(B).", call. = FALSE)
    w <- as.numeric(weights)
  }

  ## scale the WEIGHTS (not the output) so the invariant f == B %*% w holds
  if (!is.null(amplitude)) {
    f0 <- as.numeric(B %*% w)
    w  <- w * (amplitude / .safe_sd(f0))
  }

  f <- as.numeric(B %*% w)

  list(basis = B, weights = w, fitted = f)
}


#' Simulate Hierarchical Smooth Effects from a Basis Matrix
#'
#' @description
#' Generates group-specific nonlinear functions that share a common basis and a
#' shared population mean curve, with group-specific deviations around it.
#'
#' The population curve follows a first-order random walk (or is supplied via
#' `weights`):
#' \deqn{\bar w_1 = \tau z_1, \quad \bar w_k = \bar w_{k-1} + \tau z_k.}
#'
#' Group coefficients are \eqn{w_{jk} = \bar w_k + \delta_{jk}}, where the
#' deviations \eqn{\delta_{jk}} are controlled by `deviation`:
#' \itemize{
#'   \item `"iid"` (default): \eqn{\delta_{jk} \sim N(0, \sigma_{group})} --
#'     independent across basis functions (each group's deviation is NOT itself
#'     smooth; the group curve is the smooth mean plus rough wiggle).
#'   \item `"rw"` (the refinement): each group's deviation is itself a
#'     first-order random walk anchored at zero,
#'     \eqn{\delta_{j1} = 0,\ \delta_{jk} = \delta_{j,k-1} + \sigma_{group} z_{jk}},
#'     so each group's curve is smooth in its own right. The anchor at 0 keeps the
#'     population curve owning the overall level (identifiability).
#' }
#'
#' If `amplitude` is not `NULL`, a SINGLE common scale factor (derived from the
#' population curve) is applied to the population weights and all group weights,
#' so the group curves still scatter correctly around the population curve. The
#' invariant `fitted[i] == B[i,] %*% group_weights[group_i,]` is preserved.
#'
#' @param B Numeric basis matrix (rows = observations, columns = basis funs).
#' @param group Vector of group memberships, one per row of `B`.
#' @param weights Optional population coefficients (length `ncol(B)`). If `NULL`,
#'   generated by a random walk.
#' @param tau Random-walk step size for the population curve.
#' @param sigma_group SD of the group deviations.
#' @param deviation Either `"iid"` (rough deviations) or `"rw"` (smooth
#'   deviations; the refinement).
#' @param amplitude Target SD for the population curve, applied as a common scale
#'   to all weights. `NULL` leaves the raw scale.
#'
#' @return List with `basis`, `population_weights`, `group_weights`
#'   (J x K matrix of TRUE group coefficients), `population_curve`, and `fitted`.
#'   `fitted` equals the row-wise `B %*% group_weights[group,]`.
#'
#' @examples
#' x <- sort(runif(400)); group <- rep(1:4, each = 100)
#' B <- splines::bs(x, df = 8)
#' f <- simulate_hierarchical_spline_effect(B, group, deviation = "rw")
#' @export
simulate_hierarchical_spline_effect <- function(
    B,
    group,
    weights     = NULL,
    tau         = 0.15,
    sigma_group = 0.5,
    deviation   = c("iid", "rw"),
    amplitude   = 5
){
  deviation <- match.arg(deviation)
  K <- ncol(B)
  groups <- sort(unique(group))
  J <- length(groups)

  ## population weights
  if (is.null(weights)) {
    w_bar <- .rw1(K, tau)
  } else {
    if (length(weights) != K)
      stop("'weights' must have length equal to ncol(B).", call. = FALSE)
    w_bar <- as.numeric(weights)
  }

  ## group weights = population + deviation
  W <- matrix(NA_real_, nrow = J, ncol = K)
  for (j in seq_len(J)) {
    if (deviation == "iid") {
      delta <- stats::rnorm(K, sd = sigma_group)         # rough, independent
    } else {                                             # "rw": smooth, anchored at 0
      delta <- numeric(K)                                 # delta[1] = 0 (anchor)
      z <- stats::rnorm(K)
      if (K > 1) for (k in 2:K) delta[k] <- delta[k - 1] + sigma_group * z[k]
    }
    W[j, ] <- w_bar + delta
  }

  ## ONE common scale (from the population curve), applied to ALL weights,
  ## so group curves keep their correct scatter around the population curve.
  if (!is.null(amplitude)) {
    pop0  <- as.numeric(B %*% w_bar)
    scale <- amplitude / .safe_sd(pop0)
    w_bar <- w_bar * scale
    W     <- W * scale
  }

  ## evaluate
  population_curve <- as.numeric(B %*% w_bar)
  fitted <- numeric(nrow(B))
  for (j in seq_len(J)) {
    idx <- group == groups[j]
    fitted[idx] <- B[idx, , drop = FALSE] %*% W[j, ]
  }

  list(
    basis              = B,
    population_weights = w_bar,
    group_weights      = W,
    population_curve   = population_curve,
    fitted             = fitted
  )
}


#' Simulate a Nonlinear Effect Using a P-Spline
#'
#' @description
#' Builds a B-spline basis from `x` and samples smooth coefficients via a
#' first-order random-walk prior -- i.e. a P-spline. The random walk supplies the
#' smoothness penalty; `bs()` supplies only the basis.
#'
#' \deqn{w_1 = \tau z_1, \quad w_k = w_{k-1} + \tau z_k, \quad z_k \sim N(0,1).}
#'
#' If `effect_size` is not `NULL`, the WEIGHTS are scaled so `sd(fitted)`
#' equals it, keeping `fitted == B %*% weights`.
#'
#' @param x Predictor values.
#' @param df Number of basis functions.
#' @param tau Random-walk step size (smaller = smoother).
#' @param effect_size Target SD of the effect (scales the weights). `NULL` to
#'   leave the raw scale.
#' @param boundary_knots Boundary knots for the basis; defaults to `range(x)`.
#'
#' @return List with `basis`, `weights` (TRUE coefficients), `fitted`.
#'
#' @examples
#' x <- seq(0, 10, length.out = 200)
#' s <- simulate_p_spline(x, df = 8, tau = 0.2)
#' @export
simulate_p_spline <- function(
    x,
    df             = 8,
    tau            = 0.2,
    effect_size    = 1,
    boundary_knots = range(x)
){
  B <- splines::bs(x, df = df, intercept = TRUE, Boundary.knots = boundary_knots)
  K <- ncol(B)

  w <- .rw1(K, tau)

  if (!is.null(effect_size)) {
    f0 <- as.numeric(B %*% w)
    w  <- w * (effect_size / .safe_sd(f0))
  }

  f <- as.numeric(B %*% w)
  list(basis = B, weights = w, fitted = f)
}


#' Simulate Group-Specific P-Splines (independent OR hierarchical)
#'
#' @description
#' Generates a smooth nonlinear function per group over a shared B-spline basis.
#'
#' \strong{Important:} the original version of this function drew a fresh,
#' independent random walk per group with \emph{no shared population curve} -- so
#' groups were NOT pooled and it was not hierarchical despite the name. This
#' version makes the behaviour explicit via `pool`:
#' \itemize{
#'   \item `pool = FALSE` (the old behaviour): each group is an INDEPENDENT smooth
#'     random walk. No population mean, no shrinkage.
#'   \item `pool = TRUE` (default): a shared population random walk `w_bar` is
#'     drawn first, and each group deviates around it by `sigma_group` (a genuine
#'     hierarchical structure).
#' }
#'
#' When `pool = TRUE`, group weights are \eqn{w_{jk} = \bar w_k + \delta_{jk}}
#' with `deviation` controlling whether the deviations are `"iid"` (rough) or
#' `"rw"` (smooth, anchored at 0), exactly as in
#' `simulate_hierarchical_spline_effect`.
#'
#' @param x Predictor values.
#' @param group Group membership per observation.
#' @param df Number of basis functions.
#' @param tau Random-walk step size.
#' @param sigma_group SD of group deviations (used only when `pool = TRUE`).
#' @param pool Logical; `TRUE` for hierarchical (shared mean + deviations),
#'   `FALSE` for independent per-group walks.
#' @param deviation `"iid"` or `"rw"` (used only when `pool = TRUE`).
#' @param effect_size Target SD applied as a common scale. `NULL` to leave raw.
#'
#' @return List with `basis`, `population_weights` (or `NULL` if unpooled),
#'   `group_weights` (J x K), and `fitted`.
#'
#' @examples
#' x <- runif(300, 0, 100); group <- rep(1:10, each = 30)
#' s <- simulate_group_p_spline(x, group, df = 8, pool = TRUE)
#' @export
simulate_group_p_spline <- function(
    x,
    group,
    df          = 8,
    tau         = 0.2,
    sigma_group = 0.5,
    pool        = TRUE,
    deviation   = c("iid", "rw"),
    effect_size = 1
){
  deviation <- match.arg(deviation)
  B <- splines::bs(x, df = df, intercept = TRUE)
  K <- ncol(B)
  groups <- sort(unique(group))
  J <- length(groups)

  W <- matrix(NA_real_, nrow = J, ncol = K)
  w_bar <- NULL

  if (pool) {
    w_bar <- .rw1(K, tau)                                  # shared population curve
    for (j in seq_len(J)) {
      if (deviation == "iid") {
        delta <- stats::rnorm(K, sd = sigma_group)
      } else {
        delta <- numeric(K)                                # anchored at 0
        z <- stats::rnorm(K)
        if (K > 1) for (k in 2:K) delta[k] <- delta[k - 1] + sigma_group * z[k]
      }
      W[j, ] <- w_bar + delta
    }
  } else {
    for (j in seq_len(J)) W[j, ] <- .rw1(K, tau)           # independent per group
  }

  ## common scale
  if (!is.null(effect_size)) {
    ref   <- if (pool) as.numeric(B %*% w_bar) else as.numeric(B %*% W[1, ])
    scale <- effect_size / .safe_sd(ref)
    W     <- W * scale
    if (pool) w_bar <- w_bar * scale
  }

  ## evaluate
  f <- numeric(length(x))
  for (j in seq_len(J)) {
    idx <- which(group == groups[j])
    f[idx] <- B[idx, , drop = FALSE] %*% W[j, ]
  }

  list(basis = B, population_weights = w_bar, group_weights = W, fitted = f)
}


#' Create Gaussian Radial Basis Functions
#'
#' @description
#' Localized Gaussian bump basis. Column k is
#' \eqn{B_k(x) = \exp(-\tfrac12((x - \kappa_k)/\ell)^2)}.
#'
#' @param x Numeric predictor.
#' @param knots Either the number of evenly-spaced centres, or a vector of centre
#'   locations.
#' @param length_scale Bump width \eqn{\ell}. If `NULL`, defaults to
#'   `overlap * (knot spacing)`, giving adjacent bumps that overlap sensibly.
#' @param overlap Multiplier on the knot spacing for the default width (default
#'   1.5, so bumps overlap rather than sit isolated).
#'
#' @return An `n x n_knots` matrix (always a matrix, even for `n = 1`), with
#'   column names and a `knots`/`length_scale` attribute for reuse when fitting.
#'
#' @examples
#' x <- runif(200, 0, 10)
#' B <- make_radial_basis(x, knots = 8)
#' @export
make_radial_basis <- function(
    x,
    knots        = 8,
    length_scale = NULL,
    overlap      = 1.5
){
  if (length(knots) == 1) {
    knots <- seq(min(x), max(x), length.out = knots)
  }
  n_knots <- length(knots)

  if (is.null(length_scale)) {
    spacing <- if (n_knots > 1) mean(diff(sort(knots))) else diff(range(x))
    length_scale <- overlap * spacing
  }

  ## build columns; force a matrix even when length(x) == 1
  B <- vapply(
    knots,
    function(kappa) exp(-0.5 * ((x - kappa) / length_scale)^2),
    numeric(length(x))
  )
  B <- matrix(B, nrow = length(x), ncol = n_knots)
  colnames(B) <- paste0("rbf_", seq_len(n_knots))

  attr(B, "knots")        <- knots
  attr(B, "length_scale") <- length_scale
  B
}