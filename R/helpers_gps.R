#### Simulation helpers — exact GP ####

#' Simulate a 1-D Gaussian Process draw
#'
#' Draws one realization of a mean-zero Gaussian Process with an
#' exponentiated quadratic (squared-exponential) kernel over a 1-D input,
#' via the exact (dense) Cholesky construction. Non-centered: an iid
#' standard normal vector `z` is transformed by the kernel's Cholesky
#' factor, `f = L %*% z`.
#'
#' @param x Numeric vector of 1-D input locations (e.g. time, dose, age).
#' @param eta Positive amplitude (marginal standard deviation) of the GP.
#' @param rho Positive length-scale of the GP, on the same scale as `x`.
#'
#' @return Numeric vector, same length as `x`: one draw of the latent GP
#'   function `f` evaluated at `x`. No observation noise is added — this
#'   is the noiseless truth surface, not the data.
#'
#' @examples
#' x <- sort(runif(50, 0, 10))
#' f <- simulate_gp_1d(x, eta = 1, rho = 2)
#' plot(x, f, type = "l")
#'
#' @export
simulate_gp_1d <- function(x, eta = 1.5, rho = 2.5) {
  n <- length(x)

  # Distance matrix
  D <- as.matrix(dist(matrix(x, ncol = 1)))

  # Kernel
  K <- eta^2 * exp(-D^2 / (2 * rho^2))
  diag(K) <- diag(K) + 1e-6

  # Cholesky and draw the latent surface
  L <- t(chol(K))
  z <- rnorm(n, 0, 1)
  f <- as.numeric(L %*% z)

  return(f)
}

#' Simulate a 2-D Gaussian Process draw over coordinates
#'
#' Draws one realization of a mean-zero Gaussian Process with an
#' exponentiated quadratic kernel over a 2-D input (e.g. longitude/
#' latitude), via the exact (dense) Cholesky construction. Same pattern as
#' simulate_gp_1d(), with Euclidean distance in place of absolute
#' difference.
#'
#' @param longitude Numeric vector of longitude (or x) coordinates.
#' @param latitude Numeric vector of latitude (or y) coordinates, same
#'   length as `longitude`.
#' @param eta Positive amplitude (marginal standard deviation) of the GP.
#' @param rho Positive length-scale of the GP, on the same scale as the
#'   coordinates.
#'
#' @return Numeric vector, one draw of the latent GP surface `f` evaluated
#'   at each coordinate pair. No observation noise is added.
#'
#' @details Dense (exact) construction: cost scales as O(n^3) via the
#'   Cholesky factorization, so this is only practical up to roughly a few
#'   thousand points. For larger `n`, use simulate_hsgp_2d().
#'
#' @examples
#' lon <- runif(80, 0, 10)
#' lat <- runif(80, 0, 10)
#' f <- simulate_gp_2d(lon, lat, eta = 1.5, rho = 2.5)
#'
#' @export
simulate_gp_2d <- function(longitude, latitude, eta = 1.5, rho = 2.5) {
  stopifnot(length(longitude) == length(latitude))
  n <- length(longitude)

  # Bind the coords
  coord <- cbind(longitude, latitude)

  # Build a distance matrix
  D <- as.matrix(dist(coord))

  # Turn the distance matrix into a kernel
  K <- eta^2 * exp(-D^2 / (2 * rho^2))
  diag(K) <- diag(K) + 1e-6

  # Cholesky and draw the latent surface
  L <- t(chol(K))
  z <- rnorm(n, 0, 1)
  f <- as.numeric(L %*% z)

  return(f)
}

#### Simulation helpers — Hilbert Space GP approximation ####

#' Simulate a 1-D Hilbert Space GP (HSGP) approximation draw
#'
#' Approximates a mean-zero GP with an exponentiated quadratic kernel using
#' a fixed Laplacian-eigenfunction basis on a bounded domain, instead of
#' the dense n x n kernel. Cost scales O(n * M) rather than O(n^3), so this
#' is the version to use once simulate_gp_1d() becomes impractical.
#'
#' @param x Numeric vector of 1-D input locations. Does not need to be
#'   pre-centered — it is centered internally so the basis domain
#'   `-L, L` is well placed around the data.
#' @param M Number of basis functions (eigenfunctions). Larger `M` gives a
#'   more accurate approximation at higher cost; check sensitivity by
#'   increasing `M` and confirming `f` stops changing appreciably.
#' @param rho Positive length-scale of the approximated GP.
#' @param alpha Positive amplitude (marginal standard deviation) of the
#'   approximated GP.
#' @param c Boundary factor: how far the basis domain extends past the
#'   range of centered `x`, as a multiple of `max(abs(x - mean(x)))`.
#'   Typical values are 1.2--2; too small biases the fit near the edges of
#'   the data.
#'
#' @return Numeric vector, same length as `x`: one draw of the latent
#'   approximate GP surface `f`. No observation noise is added.
#'
#' @seealso simulate_gp_1d() for the exact (dense) equivalent, useful for
#'   validating the approximation on the same synthetic input.
#'
#' @examples
#' x <- sort(runif(500, 0, 10))
#' f <- simulate_hsgp_1d(x, M = 20, rho = 2, alpha = 1)
#'
#' @export
simulate_hsgp_1d <- function(x, M = 15, rho = 2, alpha = 1, c = 1.5) {
  n <- length(x)

  # Center so the basis domain [-L, L] is placed around the data
  x_c <- x - mean(x)

  ## Computational domain
  L <- c * max(abs(x_c))

  ## 1-D basis (Laplacian eigenfunctions on [-L, L])
  Phi <- sapply(1:M, function(j) {
    sin(j * pi * (x_c + L) / (2 * L)) / sqrt(L)
  })

  ## Eigenvalues
  lambda <- (seq_len(M) * pi / (2 * L))^2

  ## Spectral density (exponentiated quadratic kernel, D = 1)
  sqrt_spd <- alpha * (2 * pi)^(1 / 4) * sqrt(rho) * exp(-0.25 * rho^2 * lambda)

  ## Basis coefficients
  z <- rnorm(M, 0, 1)

  ## Latent GP
  f <- as.vector(Phi %*% (sqrt_spd * z))

  return(f)
}

#' Simulate a 2-D Hilbert Space GP (HSGP) approximation draw
#'
#' Approximates a mean-zero GP over 2-D coordinates with an exponentiated
#' quadratic kernel using a tensor product of 1-D Laplacian-eigenfunction
#' bases, instead of the dense n x n kernel. Cost scales roughly
#' O(n * M1 * M2) rather than O(n^3), so this is the version to use once
#' simulate_gp_2d() becomes impractical (n in the thousands).
#'
#' @param longitude Numeric vector of longitude (or x) coordinates.
#' @param latitude Numeric vector of latitude (or y) coordinates, same
#'   length as `longitude`. Neither needs to be pre-centered — both are
#'   centered internally.
#' @param M1,M2 Number of basis functions along the longitude and latitude
#'   axes respectively. Total basis size is `M1 * M2`; check sensitivity by
#'   increasing both and confirming `f` stops changing appreciably.
#' @param rho Positive length-scale, shared isotropically across both axes.
#' @param alpha Positive amplitude (marginal standard deviation) of the
#'   approximated GP.
#' @param c Boundary factor, same role as in simulate_hsgp_1d(), applied
#'   independently on each axis.
#'
#' @return Numeric vector, one draw of the latent approximate GP surface
#'   `f` evaluated at each coordinate pair. No observation noise is added
#'   — add it yourself (e.g. `y <- f + rnorm(n, 0, sigma)`) so the same
#'   noise-free `f` can be compared directly against simulate_gp_2d()
#'   when validating the approximation.
#'
#' @seealso simulate_gp_2d() for the exact (dense) equivalent.
#'
#' @examples
#' lon <- runif(2000, 0, 10)
#' lat <- runif(2000, 0, 10)
#' f <- simulate_hsgp_2d(lon, lat, M1 = 15, M2 = 15, rho = 2, alpha = 1)
#'
#' @export
simulate_hsgp_2d <- function(longitude, latitude, M1 = 15, M2 = 15, rho = 2, alpha = 1, c = 1.5) {
  stopifnot(length(longitude) == length(latitude))
  n <- length(longitude)

  # Center so each basis domain [-L, L] is placed around the data
  lon_c <- longitude - mean(longitude)
  lat_c <- latitude - mean(latitude)

  ## Computational domain (independent per axis)
  Lx <- c * max(abs(lon_c))
  Ly <- c * max(abs(lat_c))

  ## 1-D bases
  Px <- sapply(1:M1, function(j) sin(j * pi * (lon_c + Lx) / (2 * Lx)) / sqrt(Lx))
  Py <- sapply(1:M2, function(j) sin(j * pi * (lat_c + Ly) / (2 * Ly)) / sqrt(Ly))

  ## Tensor-product basis
  M      <- M1 * M2
  Phi    <- matrix(NA, n, M)
  lambda <- numeric(M)
  m <- 1
  for (j1 in 1:M1) {
    for (j2 in 1:M2) {
      Phi[, m]  <- Px[, j1] * Py[, j2]
      lambda[m] <- (j1 * pi / (2 * Lx))^2 + (j2 * pi / (2 * Ly))^2
      m <- m + 1
    }
  }

  ## Spectral density (exponentiated quadratic kernel, D = 2, isotropic rho)
  sqrt_spd <- alpha * sqrt(2 * pi) * rho * exp(-0.25 * rho^2 * lambda)

  ## Basis coefficients
  z <- rnorm(M, 0, 1)

  ## Latent GP
  f <- as.vector(Phi %*% (sqrt_spd * z))

  return(f)
}