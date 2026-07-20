
#' Simulate independent latent variable(s)
#'
#' Generative process (Statistical Rethinking notation):
#'   eta_i ~ Normal(mu, sigma)     for i = 1..n
#'
#' Non-centered form (default):
#'   z_i   ~ Normal(0, 1)
#'   eta_i  = mu + sigma * z_i
#'
#' @param n Number of units (rows).
#' @param mean Mean mu of the latent distribution.
#' @param sd Standard deviation sigma of the latent distribution.
#' @param ncp Logical; TRUE (default) draws via the non-centered
#'   parameterization eta_i = mu + sigma * z_i.
#'
#' @return list(
#'   data  = list(eta = numeric vector, length n),
#'   truth = list(mean = mean, sd = sd, n = n)
#' )
#' @examples
#' sim <- simulate_latent(n = 500, mean = 0, sd = 1)
#' str(sim$data$eta)
#' @export
simulate_latent <- function(n, mean = 0, sd = 1, ncp = TRUE) {
  stopifnot(n > 0, sd > 0)

  if (ncp) {
    z   <- rnorm(n, 0, 1)
    eta <- mean + sd * z
  } else {
    eta <- rnorm(n, mean, sd)
  }

  list(
    data  = list(eta = eta),
    truth = list(mean = mean, sd = sd, n = n)
  )
}


#' Simulate correlated latent variables
#'
#' Generative process:
#'   (eta_1, ..., eta_K)_i ~ MVN(mu, Sigma)     for i = 1..n
#'
#' Internally uses the Cholesky factor L of Sigma (Sigma = L L^T) and
#' draws via the multivariate non-centered parameterization:
#'   z_i    ~ Normal(0, I_K)
#'   Eta_i   = mu + L %*% z_i
#' This is the multivariate analogue of simulate_latent()'s NCP and
#' avoids a hard dependency on MASS::mvrnorm.
#'
#' Supply either `Sigma` directly, or `sds` + `Rho` (spread/wobble +
#' correlation), which is usually the more interpretable way to reason
#' about correlated constructs.
#'
#' @param n Number of units (rows).
#' @param Sigma K x K covariance matrix. Either Sigma or (sds + Rho)
#'   must be supplied.
#' @param sds Optional length-K vector of standard deviations, used
#'   with Rho to build Sigma = diag(sds) %*% Rho %*% diag(sds).
#' @param Rho Optional K x K correlation matrix, paired with sds.
#' @param means Length-K vector of latent means (default: zeros).
#' @param names Optional length-K character vector naming the latents
#'   (default: eta1, eta2, ...).
#'
#' @return list(
#'   data  = list(Eta = n x K matrix of latent draws),
#'   truth = list(means = ..., Sigma = ..., Rho = ..., sds = ..., L = ..., K = ..., n = ...)
#' )
#' @examples
#' Rho <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
#' sim <- simulate_correlated_latents(n = 500, sds = c(1, 1), Rho = Rho,
#'                                     names = c("vigour", "stress_tolerance"))
#' @export
simulate_correlated_latents <- function(n, Sigma = NULL, sds = NULL, Rho = NULL,
                                         means = NULL, names = NULL) {
  stopifnot(n > 0)

  if (is.null(Sigma)) {
    stopifnot(!is.null(sds), !is.null(Rho))
    K <- length(sds)
    Sigma <- diag(sds, K, K) %*% Rho %*% diag(sds, K, K)
  } else {
    K <- nrow(Sigma)
    if (is.null(Rho)) {
      sds <- sqrt(diag(Sigma))
      Rho <- diag(1 / sds, K, K) %*% Sigma %*% diag(1 / sds, K, K)
    }
  }

  if (is.null(means)) means <- rep(0, K)
  if (is.null(names)) names <- paste0("eta", seq_len(K))

  L <- t(chol(Sigma))                              # Sigma = L L^T
  Z <- matrix(rnorm(n * K), nrow = K, ncol = n)     # z_i ~ Normal(0, I_K)
  Eta <- t(means + L %*% Z)                         # n x K
  colnames(Eta) <- names

  list(
    data  = list(Eta = Eta),
    truth = list(means = means, Sigma = Sigma, Rho = Rho,
                 sds = sqrt(diag(Sigma)), L = L, K = K, n = n)
  )
}


#' Simulate a hierarchical (multilevel) latent variable
#'
#' Generative process:
#'   mu_j   ~ Normal(0, sigma_mu)          j = 1..J   (group means)
#'   eta_i  ~ Normal(mu_J i, sigma)       i = 1..n   (unit latent)
#'
#' Non-centered parameterization (default):
#'   mu_j   = sigma_mu * z_mu_j,   z_mu_j ~ Normal(0, 1)
#'   eta_i  = mu_J i + sigma * z_i,  z_i ~ Normal(0, 1)
#'
#' Example structure: Hospital -> Patient latent quality.
#'
#' @param n_groups Number of groups J.
#' @param n_per_group Single integer (equal group sizes) or an integer
#'   vector of length n_groups.
#' @param sigma_mu SD of the group means (spread across groups).
#' @param sigma Residual SD of units within a group (wobble).
#' @param ncp Logical; use the non-centered parameterization (default TRUE).
#'
#' @return list(
#'   data  = list(eta = numeric vector, group_id = integer vector),
#'   truth = list(mu_j = ..., sigma_mu = ..., sigma = ..., n_groups = ..., n_per_group = ...)
#' )
#' @examples
#' sim <- simulate_hierarchical_latent(n_groups = 20, n_per_group = 15,
#'                                      sigma_mu = 1, sigma = 0.5)
#' @export
simulate_hierarchical_latent <- function(n_groups, n_per_group,
                                          sigma_mu = 1, sigma = 1, ncp = TRUE) {
  stopifnot(n_groups > 0, sigma_mu > 0, sigma > 0)

  if (length(n_per_group) == 1) {
    n_per_group <- rep(n_per_group, n_groups)
  }
  stopifnot(length(n_per_group) == n_groups)

  group_id <- rep(seq_len(n_groups), times = n_per_group)
  n <- length(group_id)

  if (ncp) {
    z_mu <- rnorm(n_groups, 0, 1)
    mu_j <- sigma_mu * z_mu
    z    <- rnorm(n, 0, 1)
    eta  <- mu_j[group_id] + sigma * z
  } else {
    mu_j <- rnorm(n_groups, 0, sigma_mu)
    eta  <- rnorm(n, mu_j[group_id], sigma)
  }

  list(
    data  = list(eta = eta, group_id = group_id),
    truth = list(mu_j = mu_j, sigma_mu = sigma_mu, sigma = sigma,
                 n_groups = n_groups, n_per_group = n_per_group)
  )
}


#' Simulate a latent variable predicted by observed covariates (MIMIC-style)
#'
#' Generative process:
#'   eta_i  = beta0 + beta_1 x_1i + beta_2 x_2i + ... + zeta_i
#'   zeta_i ~ Normal(0, sigma_zeta)
#'
#' Example structure: Training, Age, Sleep -> Fitness (latent).
#'
#' @param X Numeric matrix or data.frame of covariates (n x P). Columns
#'   should already be centered/scaled as desired; this function does
#'   not transform them.
#' @param betas Numeric vector of length P, one coefficient per column of X.
#' @param beta0 Intercept.
#' @param sigma_zeta Residual SD of the structural disturbance zeta.
#'
#' @return list(
#'   data  = list(eta = numeric vector, X = matrix),
#'   truth = list(beta0 = ..., betas = ..., sigma_zeta = ..., n = ...)
#' )
#' @examples
#' X   <- data.frame(age = rnorm(300), sleep = rnorm(300))
#' sim <- simulate_latent_regression(X, betas = c(0.5, -0.3), beta0 = 0, sigma_zeta = 1)
#' @export
simulate_latent_regression <- function(X, betas, beta0 = 0, sigma_zeta = 1) {
  X <- as.matrix(X)
  stopifnot(ncol(X) == length(betas), sigma_zeta > 0)

  n    <- nrow(X)
  zeta <- rnorm(n, 0, sigma_zeta)
  eta  <- as.numeric(beta0 + X %*% betas + zeta)

  list(
    data  = list(eta = eta, X = X),
    truth = list(beta0 = beta0, betas = betas, sigma_zeta = sigma_zeta, n = n)
  )
}


#' Simulate one latent variable as a structural path from another
#'
#' Generative process:
#'   eta2_i = beta * eta1_i + zeta_i
#'   zeta_i ~ Normal(0, sigma_zeta)
#'
#' Example structure: Teaching quality -> Student engagement.
#'
#' @param eta1 Numeric vector, the upstream (predictor) latent variable.
#' @param beta Structural path coefficient eta1 -> eta2.
#' @param sigma_zeta Residual SD of the structural disturbance zeta.
#'
#' @return list(
#'   data  = list(eta1 = eta1, eta2 = numeric vector),
#'   truth = list(beta = beta, sigma_zeta = sigma_zeta, n = n)
#' )
#' @examples
#' teaching <- simulate_latent(n = 400)$data$eta
#' sim      <- simulate_latent_path(teaching, beta = 0.6, sigma_zeta = 0.8)
#' @export
simulate_latent_path <- function(eta1, beta, sigma_zeta = 1) {
  stopifnot(sigma_zeta > 0)

  n    <- length(eta1)
  zeta <- rnorm(n, 0, sigma_zeta)
  eta2 <- beta * eta1 + zeta

  list(
    data  = list(eta1 = eta1, eta2 = eta2),
    truth = list(beta = beta, sigma_zeta = sigma_zeta, n = n)
  )
}