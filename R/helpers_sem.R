
#' Simulate an outcome from a latent (or observed) predictor
#'
#' Generative process
#'   mu_i = beta0 + beta1 * eta_i  
#'   gaussian:  Y_i ~ Normal(mu_i, sigma)
#'   bernoulli: Y_i ~ Bernoulli(logistic(mu_i))
#'   poisson:   Y_i ~ Poisson(exp(mu_i))
#'
#' The latent enters through the linear predictor and then a
#' family-specific inverse link and sampling distribution, e.g.
#' Engagement -> Exam score, or Disease severity -> Mortality.
#'
#' @param eta Numeric vector, the latent (or observed) predictor.
#' @param family One of "gaussian", "bernoulli", "poisson".
#' @param beta0 Intercept.
#' @param beta1 Coefficient on eta.
#' @param sigma Residual SD, required only when family == "gaussian".
#' @param X Optional additional covariate matrix (n x P) entering the
#'   linear predictor additively (for outcome models with more than
#'   one predictor).
#' @param betas_x Optional numeric vector of length P, coefficients for X.
#'
#' @return list(
#'   data  = list(Y = numeric/integer vector, eta = eta, X = X),
#'   truth = list(beta0 = ..., beta1 = ..., sigma = ..., betas_x = ...,
#'                 family = ..., n = ...)
#' )
#' @examples
#' engagement <- simulate_latent(n = 400)$data$eta
#' sim <- simulate_outcome(engagement, family = "gaussian",
#'                          beta0 = 50, beta1 = 5, sigma = 8)
#'
#' severity <- simulate_latent(n = 400)$data$eta
#' sim2 <- simulate_outcome(severity, family = "bernoulli",
#'                           beta0 = -2, beta1 = 1.5)
#' @export
simulate_outcome <- function(eta, family = c("gaussian", "bernoulli", "poisson"),
                              beta0, beta1, sigma = NULL, X = NULL, betas_x = NULL) {
  family <- match.arg(family)
  n <- length(eta)

  mu <- beta0 + beta1 * eta
  if (!is.null(X)) {
    stopifnot(!is.null(betas_x))
    X  <- as.matrix(X)
    mu <- mu + as.numeric(X %*% betas_x)
  }

  Y <- switch(family,
    gaussian = {
      stopifnot(!is.null(sigma), sigma > 0)
      rnorm(n, mu, sigma)
    },
    bernoulli = {
      pi_hat <- plogis(mu)
      rbinom(n, 1, pi_hat)
    },
    poisson = {
      lambda_hat <- exp(mu)
      rpois(n, lambda_hat)
    }
  )

  list(
    data  = list(Y = Y, eta = eta, X = X),
    truth = list(beta0 = beta0, beta1 = beta1, sigma = sigma,
                 betas_x = betas_x, family = family, n = n)
  )
}


#' Check basic identification conditions for a CFA / SEM measurement model
#'
#' This function does not simulate any data. It runs a small set of
#' sanity checks on the *scale-setting* and *counting* rules used to
#' identify a latent variable's measurement model, and reports
#' warnings when a rule looks violated. Meant as a pre-flight check
#' before writing the corresponding Stan/ulam model.
#'
#' Scale-setting rule: exactly ONE of the following must hold for each
#' latent variable -- never both, never neither:
#'   (a) Var(eta) = 1                     (latent variance fixed)
#'   (b) lambda_1 = 1 for one marker item (marker-loading fixed)
#'
#' Counting rule (rule of thumb): a single latent variable needs at
#' least 3 indicators to be identified without extra equality
#' constraints. With exactly 2 indicators it is identified only under
#' extra constraints (e.g. equal loadings and/or equal error
#' variances), which this function flags but cannot verify on its own.
#' With 0 or 1 indicators it cannot be identified at all.
#'
#' @param n_indicators Integer, or an integer vector (one entry per
#'   latent variable, optionally named) giving how many indicators
#'   load on each latent.
#' @param latent_variance_fixed Logical, or a vector matching
#'   n_indicators; TRUE if Var(eta) has been fixed (typically to 1)
#'   for that latent.
#' @param marker_loading_fixed Logical, or a vector matching
#'   n_indicators; TRUE if one loading has been fixed (typically to 1)
#'   for that latent.
#' @param verbose Logical; print human-readable messages (default TRUE).
#'
#' @return Invisibly, a data.frame with one row per latent variable:
#'   latent, n_indicators, latent_variance_fixed, marker_loading_fixed,
#'   status ("ok", "underidentified", "overconstrained",
#'   "needs_constraints"), and message.
#' @examples
#' check_identification(n_indicators = 3, latent_variance_fixed = TRUE,
#'                       marker_loading_fixed = FALSE)
#'
#' check_identification(n_indicators = c(vigour = 3, stress = 2),
#'                       latent_variance_fixed = c(FALSE, TRUE),
#'                       marker_loading_fixed = c(TRUE, TRUE))
#' @export
check_identification <- function(n_indicators, latent_variance_fixed = FALSE,
                                  marker_loading_fixed = FALSE, verbose = TRUE) {

  K <- length(n_indicators)
  latent_variance_fixed <- rep(latent_variance_fixed, length.out = K)
  marker_loading_fixed  <- rep(marker_loading_fixed,  length.out = K)

  lat_names <- names(n_indicators)
  if (is.null(lat_names)) lat_names <- paste0("eta", seq_len(K))

  results <- data.frame(
    latent = lat_names,
    n_indicators = as.integer(n_indicators),
    latent_variance_fixed = latent_variance_fixed,
    marker_loading_fixed = marker_loading_fixed,
    status = NA_character_,
    message = NA_character_,
    stringsAsFactors = FALSE
  )

  for (k in seq_len(K)) {
    J    <- results$n_indicators[k]
    lvf  <- latent_variance_fixed[k]
    mlf  <- marker_loading_fixed[k]
    msgs <- character(0)
    status <- "ok"

    # --- scale-setting rule ---
    if (lvf && mlf) {
      status <- "overconstrained"
      msgs <- c(msgs, "Both latent variance AND a marker loading are fixed -- scale is set twice; drop one constraint.")
    } else if (!lvf && !mlf) {
      status <- "underidentified"
      msgs <- c(msgs, "No scale-setting constraint found -- fix either Var(eta) = 1 or one marker loading = 1.")
    }

    # --- counting rule ---
    if (J < 1) {
      status <- "underidentified"
      msgs <- c(msgs, "A latent variable with 0 indicators cannot be identified.")
    } else if (J == 1) {
      status <- "underidentified"
      msgs <- c(msgs, "Only 1 indicator: loading and error variance cannot be separated without external constraints.")
    } else if (J == 2) {
      if (status == "ok") status <- "needs_constraints"
      msgs <- c(msgs, "Only 2 indicators: identified only with extra equality constraints (e.g. equal loadings and/or equal error variances); verify manually.")
    }

    if (length(msgs) == 0) msgs <- "Identification looks fine (scale set once; >= 3 indicators)."

    results$status[k]  <- status
    results$message[k] <- paste(msgs, collapse = " ")
  }

  if (verbose) {
    for (k in seq_len(K)) {
      cat(sprintf("[%s] status: %s\n  %s\n\n",
                  results$latent[k], results$status[k], results$message[k]))
    }
  }

  invisible(results)
}