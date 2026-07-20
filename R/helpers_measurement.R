
#' Simulate observed indicators from a latent variable (measurement equation)
#'
#' Generative process:
#'   y_ij = nu_j + lambda_j * eta_i + epsilon_ij
#'   epsilon_ij ~ Normal(0, sigma_j)
#'
#' Supports a single latent variable (eta as a numeric vector, loadings
#' as a vector of length J) or several latents loading onto a shared
#' item set (eta as an n x K matrix, loadings as a K x J matrix -- this
#' allows cross-loadings for multi-factor CFA).
#'
#' @param eta Numeric vector (single latent, length n) or n x K matrix
#'   (K latents), as produced by simulate_latent() / simulate_correlated_latents().
#' @param loadings Numeric vector of length J (single latent) giving
#'   lambda_j, or a K x J matrix (multiple latents) giving lambda_kj.
#' @param intercepts Numeric vector of length J, item intercepts nu_j
#'   (default: zeros).
#' @param sigma Numeric scalar or vector of length J, residual SD(s)
#'   epsilon_j (default: 1 for every item).
#' @param names Optional character vector of length J naming the items
#'   (default: item1, item2, ...).
#'
#' @return list(
#'   data  = list(Y = n x J matrix of observed indicators),
#'   truth = list(intercepts = ..., loadings = ..., sigma = ..., n = ..., J = ...)
#' )
#' @examples
#' eta <- simulate_latent(n = 500)$data$eta
#' sim <- simulate_indicators(eta, loadings = c(1, 0.8, 1.2))
#'
#' # two correlated factors, each with its own items (block-diagonal loadings)
#' Eta <- simulate_correlated_latents(n = 500, sds = c(1, 1),
#'                                     Rho = matrix(c(1, .5, .5, 1), 2))$data$Eta
#' Lambda <- matrix(c(1, 0.9, 0.8, 0,   0,   0,
#'                     0, 0,   0,   1, 0.7, 1.1), nrow = 2, byrow = TRUE)
#' sim2 <- simulate_indicators(Eta, loadings = Lambda)
#' @export
simulate_indicators <- function(eta, loadings, intercepts = NULL,
                                 sigma = NULL, names = NULL) {

  if (is.null(dim(eta))) {
    Eta    <- matrix(eta, ncol = 1)
    Lambda <- matrix(loadings, nrow = 1)
  } else {
    Eta    <- as.matrix(eta)
    Lambda <- as.matrix(loadings)
    stopifnot(nrow(Lambda) == ncol(Eta))
  }

  n <- nrow(Eta)
  J <- ncol(Lambda)

  if (is.null(intercepts)) intercepts <- rep(0, J)
  if (is.null(sigma))      sigma      <- rep(1, J)
  if (length(sigma) == 1)  sigma      <- rep(sigma, J)
  if (is.null(names))      names      <- paste0("item", seq_len(J))

  stopifnot(length(intercepts) == J, length(sigma) == J)

  lin_pred <- Eta %*% Lambda                          # n x J
  Eps <- matrix(rnorm(n * J), nrow = n, ncol = J)
  Eps <- sweep(Eps, 2, sigma, "*")

  Y <- sweep(lin_pred, 2, intercepts, "+") + Eps
  colnames(Y) <- names

  list(
    data  = list(Y = Y),
    truth = list(intercepts = intercepts, loadings = Lambda, sigma = sigma,
                 n = n, J = J)
  )
}


#' Simulate ordinal (Likert-type) indicators from a latent variable
#'
#' Generative process (ordered logit / probit measurement model)
#'   ystar_ij   = lambda_j * eta_i + epsilon_ij
#'   epsilon_ij ~ Logistic(0, 1)     ordered logit, default
#'              ~ Normal(0, 1)       ordered probit
#'   y_ij = k   iff   kappa_{j,k-1} < ystar_ij <= kappa_{j,k}
#'
#' Thresholds (kappa) cut the continuous response into K ordered
#' categories -- the same cut-point logic used elsewhere for the
#' ordered-logit family.
#'
#' @param eta Numeric vector, the latent trait (length n).
#' @param loadings Numeric vector of length J, one lambda_j per item.
#' @param thresholds Numeric vector of cut points shared across items,
#'   OR a list of length J of numeric vectors (per-item cut points).
#'   Cut points must be sorted ascending; number of categories K =
#'   length(thresholds) + 1.
#' @param dist Error distribution for the underlying continuous
#'   response: "logistic" (default, ordered logit) or "normal"
#'   (ordered probit).
#' @param names Optional character vector of length J naming the items.
#'
#' @return list(
#'   data  = list(Y = n x J integer matrix of ordinal responses,
#'                Ystar = n x J matrix of underlying continuous responses),
#'   truth = list(loadings = ..., thresholds = ..., dist = ..., n = ..., J = ...)
#' )
#' @examples
#' eta <- simulate_latent(n = 400)$data$eta
#' sim <- simulate_ordinal_indicators(eta, loadings = c(1, 1.2),
#'                                     thresholds = c(-1, 0, 1, 2))
#' table(sim$data$Y[, 1])
#' @export
simulate_ordinal_indicators <- function(eta, loadings, thresholds,
                                         dist = c("logistic", "normal"),
                                         names = NULL) {
  dist <- match.arg(dist)
  n <- length(eta)
  J <- length(loadings)

  if (!is.list(thresholds)) {
    thresholds <- rep(list(sort(thresholds)), J)
  }
  stopifnot(length(thresholds) == J)
  if (is.null(names)) names <- paste0("item", seq_len(J))

  Ystar <- matrix(NA_real_, n, J)
  Y     <- matrix(NA_integer_, n, J)

  for (j in seq_len(J)) {
    eps     <- if (dist == "logistic") rlogis(n, 0, 1) else rnorm(n, 0, 1)
    ystar_j <- loadings[j] * eta + eps
    Ystar[, j] <- ystar_j
    Y[, j] <- as.integer(cut(ystar_j,
                              breaks = c(-Inf, thresholds[[j]], Inf),
                              labels = FALSE))
  }
  colnames(Ystar) <- names
  colnames(Y)     <- names

  list(
    data  = list(Y = Y, Ystar = Ystar),
    truth = list(loadings = loadings, thresholds = thresholds, dist = dist,
                 n = n, J = J)
  )
}