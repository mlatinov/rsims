#' Simulate a Smooth Effect from a Basis Matrix
#'
#' @description
#' Generates a smooth nonlinear function by combining a basis matrix with
#' random-walk spline weights.
#'
#' This helper is independent of the type of basis used. The basis matrix may
#' come from B-splines, P-splines, radial basis functions, polynomial bases,
#' truncated power bases, or any other spline construction.
#'
#' The spline coefficients follow a first-order random walk
#'
#' \deqn{
#' w_1 = \tau z_1
#' }
#'
#' \deqn{
#' w_k = w_{k-1} + \tau z_k,
#' \qquad
#' z_k \sim Normal(0,1).
#' }
#'
#' After constructing the spline
#'
#' \deqn{
#' f(x)=Bw,
#' }
#'
#' the nonlinear effect is standardized and rescaled to a desired amplitude.
#'
#' This helper is useful when simulating nonlinear relationships because the
#' same weight-generation mechanism can be combined with many different basis
#' functions.
#'
#' @param B Numeric basis matrix. Rows correspond to observations and columns
#' correspond to basis functions.
#'
#' @param tau Positive random-walk step size controlling smoothness.
#' Smaller values produce smoother functions.
#'
#' @param amplitude Desired standard deviation of the simulated nonlinear
#' effect after rescaling.
#'
#' @return
#' A list containing
#'
#' * `basis` — the supplied basis matrix.
#' * `weights` — the simulated spline coefficients.
#' * `fitted` — the nonlinear effect evaluated for every observation.
#'
#' @examples
#' x <- sort(runif(200))
#'
#' B <- splines::bs(x, df = 8)
#'
#' f <- simulate_spline_effect(B)
#'
#' plot(x, f$fitted, type = "l")
#'
#' @export
simulate_spline_effect <- function(
    B,
    tau = 0.10,
    amplitude = 5
){

  K <- ncol(B)
  z <- rnorm(K)
  w <- numeric(K)
  w[1] <- tau * z[1]

  if(K > 1){
    for(i in 2:K){
      w[i] <- w[i - 1] + tau * z[i]
    }
  }

  f <- as.numeric(B %*% w)
  f <- f / sd(f) * amplitude

  list(
    basis   = B,
    weights = w,
    fitted  = f
  )
}
#' Simulate a Nonlinear Effect Using a P-Spline
#'
#' @description
#' Generates a smooth nonlinear function from a P-spline representation.
#' The function first constructs a B-spline basis matrix and then samples
#' spline weights using a first-order random walk prior.
#'
#' The generated nonlinear function is:
#'
#' \deqn{
#' f(x_i)=\sum_{k=1}^{K} B_k(x_i)w_k
#' }
#'
#' where \eqn{B_k(x_i)} are fixed spline basis functions and \eqn{w_k}
#' are smoothness-controlled spline weights.
#'
#' The weights are generated using a non-centered random walk:
#'
#' \deqn{
#' w_1=\tau z_1
#' }
#'
#' \deqn{
#' w_k=w_{k-1}+\tau z_k
#' }
#'
#' where:
#'
#' \itemize{
#'   \item \eqn{z_k \sim Normal(0,1)} are independent random innovations,
#'   \item \eqn{\tau} controls the smoothness of the generated curve.
#' }
#'
#' Smaller values of \eqn{\tau} create smoother functions, while larger values
#' generate more local variation.
#'
#' The resulting function is rescaled so that its standard deviation matches
#' the requested effect size.
#'
#' @param x Numeric vector containing the predictor values where the spline
#' basis is evaluated.
#'
#' @param df Number of spline basis functions.
#'
#' @param tau Random walk step size controlling smoothness.
#' Smaller values produce smoother nonlinear effects.
#'
#' @param effect_size Standard deviation of the generated nonlinear effect
#' after rescaling.
#'
#' @param boundary_knots Boundary points used for the spline basis construction.
#' Defaults to the observed range of \code{x}.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{basis}{The generated B-spline basis matrix \eqn{B}.}
#'   \item{weights}{The spline coefficients \eqn{w}.}
#'   \item{fitted}{The evaluated nonlinear function \eqn{f(x)}.}
#' }
#'
#' @examples
#' x <- seq(0, 10, length.out = 200)
#'
#' spline <- simulate_p_spline(
#'   x,
#'   df = 8,
#'   tau = 0.2
#' )
#'
#' plot(x, spline$fitted, type = "l")
#'
#' @export
simulate_p_spline <- function(
    x,
    df = 8,
    tau = 0.2,
    effect_size = 1,
    boundary_knots = range(x)
){

  # Create P-spline basis matrix
  B <- splines::bs(
    x,
    df = df,
    intercept = TRUE,
    Boundary.knots = boundary_knots
  )

  K <- ncol(B)

  # Random walk spline prior
  z <- rnorm(K)
  w <- numeric(K)
  w[1] <- tau * z[1]
  for(i in 2:K){
    w[i] <- w[i-1] + tau*z[i]
  }

  # Nonlinear function
  f <- as.vector(B %*% w)

  # Control effect magnitude
  f <- f / sd(f) * effect_size

  # Return a list with the fitted values and the parameters 
  list(
    basis = B,
    weights = w,
    fitted = f
  )
}

#' Simulate Group-Specific Random Effect P-Splines
#'
#' @description
#' Generates nonlinear functions that vary across groups using a hierarchical
#' P-spline structure.
#'
#' Each group receives its own spline coefficients, allowing each group to have
#' a different smooth relationship between the predictor and the outcome.
#'
#' The model structure is:
#'
#' \deqn{
#' f_j(x_i)=\sum_{k=1}^{K}B_k(x_i)w_{jk}
#' }
#'
#' where \eqn{j} indexes groups and \eqn{k} indexes spline basis functions.
#'
#' The group-specific spline weights follow a random walk:
#'
#' \deqn{
#' w_{j1}=\tau z_{j1}
#' }
#'
#' \deqn{
#' w_{jk}=w_{j,k-1}+\tau z_{jk}
#' }
#'
#' with:
#'
#' \deqn{
#' z_{jk}\sim Normal(0,1)
#' }
#'
#' This creates partially pooled-like nonlinear effects where all groups share
#' the same spline basis but have their own smooth deviations.
#'
#' The resulting functions can be used to simulate varying nonlinear effects
#' in hierarchical models.
#'
#' @param x Numeric vector containing predictor values.
#'
#' @param group Vector defining the group membership for each observation.
#'
#' @param df Number of spline basis functions.
#'
#' @param tau Random walk step size controlling smoothness of group curves.
#'
#' @param effect_size Standard deviation of the nonlinear effect after scaling.
#'
#' @return
#' A list containing:
#'
#' \describe{
#'   \item{basis}{The shared spline basis matrix \eqn{B}.}
#'   \item{weights}{Matrix containing group-specific spline weights.}
#'   \item{fitted}{The generated nonlinear effect for each observation.}
#' }
#'
#' @examples
#'
#' x <- runif(300, 0, 100)
#' group <- rep(1:10, each = 30)
#'
#' spline <- simulate_random_p_spline(
#'   x,
#'   group,
#'   df = 8
#' )
#'
#' plot(x, spline$fitted)
#'
#' @export
simulate_random_p_spline <- function(
    x,
    group,
    df = 8,
    tau = 0.2,
    effect_size = 1
){

  # Basis matrix shared across groups
  B <- splines::bs(x, df = df, intercept = TRUE)
  K        <- ncol(B)
  groups   <- unique(group)
  n_groups <- length(groups)

  # Storage for group spline values
  f <- numeric(length(x))

  # Storage for weights
  W <- matrix(0, nrow = n_groups, ncol = K)

  for(j in seq_along(groups)){

    # Random walk weights
    z <- rnorm(K)
    w <- numeric(K)
    w[1] <- tau*z[1]

    for(k in 2:K){
      w[k] <- w[k-1]+tau*z[k]
    }

    # save weights
    W[j,] <- w

    # evaluate spline for this group
    idx <- which(group == groups[j])
    f[idx] <- B[idx,] %*% w

  }
  # standardize nonlinear effect
  f <- f / sd(f) * effect_size

  # Return a list with fitted values and paramters 
  list(
    basis   = B,
    weights = W,
    fitted  = f
  )
}

#' Create Gaussian Radial Basis Functions
#'
#' Creates localized Gaussian bump basis functions.
#'
#' Each basis column is:
#'
#' \deqn{
#' B_k(x)=exp(-1/2((x-\kappa_k)/l)^2)
#' }
#'
#' where \eqn{\kappa_k} are knot locations.
#'
#' @param x Numeric predictor.
#' @param knots Number of radial basis centres or supplied knot values.
#' @param length_scale Width of Gaussian bumps.
#'
#' @return Matrix containing radial basis columns.
#'
#' @export
make_radial_basis <- function(
    x,
    knots = 8,
    length_scale = NULL
){

  if(length(knots)==1){
    knots <- seq(
      min(x),
      max(x),
      length.out = knots
    )
  }

  if(is.null(length_scale)){
    length_scale <- diff(range(x))/length(knots)
  }

  B <- sapply(
    knots,
    function(kappa){
      exp(-0.5*((x-kappa)/length_scale)^2)
    }
  )
  colnames(B) <- paste0("rbf_",seq_along(knots))

  return(B)
}