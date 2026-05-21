
# ---------------------------------------------------------------------------
# Matern 5/2 kernel and GP inference
#
# Main interface:
#   gp_smooth(U, tt, sigma2_n = NULL) -> list(U_star, dU_dt, hyperparams)
#
# U_star  : smoothed state  (drop-in for wendy_erts$U_star)
# dU_dt   : posterior mean of df/dt at tt  (free derivative estimate)
#
# Each state dimension is treated as an independent GP.
# Hyperparameters (sigma2, ell, sigma2_n) are optimised per dimension via
# log marginal likelihood unless supplied by the caller.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Core kernel functions
# ---------------------------------------------------------------------------

#' Matern 5/2 covariance matrix
#'
#' @param tt1,tt2 Numeric vectors of time points.
#' @param sigma2  Signal variance (> 0).
#' @param ell     Length scale (> 0).
#' @return Matrix of size length(tt1) x length(tt2).
matern52_matrix <- function(tt1, tt2, sigma2, ell) {
  tau <- outer(tt1, tt2, "-")
  r   <- abs(tau) / ell
  sigma2 * (1 + sqrt(5) * r + 5 * r^2 / 3) * exp(-sqrt(5) * r)
}

#' Matern 5/2 cross-covariance cov(f(tt1), f'(tt2))
#'
#' Derivative is with respect to the second argument.
#' @param tt1,tt2 Numeric vectors of time points.
#' @param sigma2  Signal variance (> 0).
#' @param ell     Length scale (> 0).
#' @return Matrix of size length(tt1) x length(tt2).
#' @keywords internal
matern52_Kfd <- function(tt1, tt2, sigma2, ell) {
  tau <- outer(tt1, tt2, "-")     # tau = t1 - t2
  r   <- abs(tau) / ell
  e   <- exp(-sqrt(5) * r)
  # d k(tau) / d t2 = -dk/dtau = 5*sigma2*tau/(3*ell^2) * (1 + sqrt(5)*r) * e
  (5 * sigma2 / (3 * ell^2)) * tau * (1 + sqrt(5) * r) * e
}

#' Matern 5/2 derivative kernel cov(f'(tt1), f'(tt2))
#'
#' @param tt1,tt2 Numeric vectors of time points.
#' @param sigma2  Signal variance (> 0).
#' @param ell     Length scale (> 0).
#' @return Matrix of size length(tt1) x length(tt2).
#' @keywords internal
matern52_Kdd <- function(tt1, tt2, sigma2, ell) {
  tau <- outer(tt1, tt2, "-")
  r   <- abs(tau) / ell
  e   <- exp(-sqrt(5) * r)
  # -d^2 k / dtau^2 = 5*sigma2/(3*ell^2) * (1 + sqrt(5)*r - 5*r^2) * e
  (5 * sigma2 / (3 * ell^2)) * (1 + sqrt(5) * r - 5 * tau^2 / ell^2) * e
}

#' Observation noise covariance (nugget)
#'
#' @param n        Number of observations.
#' @param sigma2_n Noise variance.
#' @return Diagonal matrix n x n.
nugget <- function(n, sigma2_n) sigma2_n * diag(n)


# ---------------------------------------------------------------------------
# Log marginal likelihood and hyperparameter optimisation
# ---------------------------------------------------------------------------

#' Log marginal likelihood for a 1-D GP with Matern 5/2
#'
#' @param log_theta Log-scale hyperparameters c(log(sigma2), log(ell), log(sigma2_n)).
#' @param tt        Observed time points (length n).
#' @param y         Observations (length n), zero-mean assumed.
#' @return Scalar log marginal likelihood.
gp_log_ml <- function(log_theta, tt, y) {
  sigma2   <- exp(log_theta[1])
  ell      <- exp(log_theta[2])
  sigma2_n <- exp(log_theta[3])
  n        <- length(tt)

  Kyy <- matern52_matrix(tt, tt, sigma2, ell) + nugget(n, sigma2_n)

  L <- tryCatch(chol(Kyy), error = function(e) NULL)
  if (is.null(L)) return(-Inf)

  alpha   <- backsolve(L, forwardsolve(t(L), y))
  log_det <- 2 * sum(log(diag(L)))

  -0.5 * (as.numeric(t(y) %*% alpha) + log_det + n * log(2 * pi))
}

#' Optimise Matern 5/2 GP hyperparameters via marginal likelihood
#'
#' Uses a grid of initial length scales (t_range times 2^(-i) for i=0..n_scales-1)
#' so that both slow and fast features are candidates, then adds random
#' restarts for robustness.
#'
#' @param tt         Observed time points.
#' @param y          Observations.
#' @param sigma2_n   Fixed noise variance; if NULL it is also optimised.
#' @param n_restarts Additional random restarts beyond the length-scale grid.
#' @param n_scales   Number of grid length scales to try (default 5).
#' @return Named list with sigma2, ell, sigma2_n.
#' @export
gp_optimize_hyperparams <- function(tt, y, sigma2_n = NULL, n_restarts = 3L,
                                    n_scales = 5L) {
  fix_noise <- !is.null(sigma2_n)
  y_sd      <- max(sd(y), 1e-6)
  t_range   <- diff(range(tt))
  n_par     <- if (fix_noise) 2L else 3L

  neg_lml <- function(log_theta) {
    lt  <- if (fix_noise) c(log_theta, log(sigma2_n)) else log_theta
    val <- gp_log_ml(lt, tt, y)
    if (!is.finite(val)) 1e10 else -val
  }

  make_init <- function(log_ell) {
    base <- c(log(y_sd^2), log_ell)
    if (!fix_noise) base <- c(base, log(y_sd^2 * 0.05)) else base
  }

  # Grid of starting length scales: t_range, t_range/2, ..., t_range/2^(n_scales-1)
  grid_lells <- log(t_range) - log(2) * seq(0, n_scales - 1L)

  inits <- lapply(grid_lells, make_init)
  for (i in seq_len(n_restarts))
    inits <- c(inits, list(make_init(log(t_range) + rnorm(1) - 1)))

  best <- list(val = Inf, par = NULL)
  for (init in inits) {
    res <- tryCatch(
      optim(init, neg_lml, method = "L-BFGS-B",
            lower = rep(-10, n_par), upper = rep(10, n_par)),
      error = function(e) list(value = Inf, par = init)
    )
    if (res$value < best$val) best <- list(val = res$value, par = res$par)
  }

  par <- best$par
  list(
    sigma2   = exp(par[1]),
    ell      = exp(par[2]),
    sigma2_n = if (fix_noise) sigma2_n else exp(par[3])
  )
}


# ---------------------------------------------------------------------------
# GP fit object and prediction
# ---------------------------------------------------------------------------

#' Fit a 1-D Matern 5/2 GP to observations
#'
#' @param tt       Training time points.
#' @param y        Observations.
#' @param sigma2   Signal variance.  NULL triggers optimisation.
#' @param ell      Length scale.     NULL triggers optimisation.
#' @param sigma2_n Noise variance.   NULL triggers optimisation.
#' @return A gp_fit object (list with kernel params and precomputed Cholesky).
#' @export
gp_fit_1d <- function(tt, y, sigma2 = NULL, ell = NULL, sigma2_n = NULL) {
  tt <- as.vector(tt)
  y  <- as.vector(y)
  n  <- length(tt)

  needs_opt <- is.null(sigma2) || is.null(ell) || is.null(sigma2_n)
  if (needs_opt) {
    hp       <- gp_optimize_hyperparams(tt, y, sigma2_n = sigma2_n)
    sigma2   <- sigma2   %||% hp$sigma2
    ell      <- ell      %||% hp$ell
    sigma2_n <- sigma2_n %||% hp$sigma2_n
  }

  Kyy <- matern52_matrix(tt, tt, sigma2, ell) + nugget(n, sigma2_n)
  L   <- chol(Kyy + 1e-10 * diag(n))   # upper triangular: Kyy = t(L) %*% L
  alpha <- backsolve(L, forwardsolve(t(L), y))

  structure(
    list(tt = tt, y = y, L = L, alpha = alpha,
         sigma2 = sigma2, ell = ell, sigma2_n = sigma2_n),
    class = "gp_fit"
  )
}

#' Predict function values from a fitted GP
#'
#' @param fit     A gp_fit object from \code{gp_fit_1d}.
#' @param tt_star Prediction time points.
#' @return List with mean (vector) and var (vector of marginal variances).
#' @export
gp_predict <- function(fit, tt_star) {
  tt_star <- as.vector(tt_star)

  Ks  <- matern52_matrix(tt_star, fit$tt, fit$sigma2, fit$ell)  # n* x n
  Kss <- matern52_matrix(tt_star, tt_star, fit$sigma2, fit$ell)

  mu  <- as.vector(Ks %*% fit$alpha)
  V   <- forwardsolve(t(fit$L), t(Ks))      # n x n*
  var <- pmax(diag(Kss) - colSums(V^2), 0)

  list(mean = mu, var = var)
}

#' Predict derivative values from a fitted GP
#'
#' Uses the analytic cross-covariance cov(f'(t*), f(t)).
#'
#' @param fit     A gp_fit object from \code{gp_fit_1d}.
#' @param tt_star Prediction time points.
#' @return List with mean (vector) and var (vector of marginal variances).
#' @export
gp_predict_deriv <- function(fit, tt_star) {
  tt_star <- as.vector(tt_star)

  # K_{df*, f}  — shape: n* x n
  Kds  <- t(matern52_Kfd(fit$tt, tt_star, fit$sigma2, fit$ell))
  # K_{df*, df*} — shape: n* x n*
  Kdds <- matern52_Kdd(tt_star, tt_star, fit$sigma2, fit$ell)

  mu  <- as.vector(Kds %*% fit$alpha)
  V   <- forwardsolve(t(fit$L), t(Kds))
  var <- pmax(diag(Kdds) - colSums(V^2), 0)

  list(mean = mu, var = var)
}


# ---------------------------------------------------------------------------
# Main WENDy interface
# ---------------------------------------------------------------------------

#' Smooth a multi-dimensional state with independent Matern 5/2 GPs
#'
#' Drop-in replacement for \code{wendy_erts}: returns \code{$U_star} (smoothed
#' state) and additionally \code{$dU_dt} (posterior mean of df/dt), which
#' \code{wendy_erts} does not provide.
#'
#' Each state dimension is fit independently. Hyperparameters are optimised
#' via marginal likelihood unless supplied via \code{hyperparams}.
#'
#' @param U           n x D matrix of noisy observations.
#' @param tt          Length-n vector of time points.
#' @param sigma2_n    Scalar or length-D vector of noise variances.
#'                    NULL triggers per-dimension optimisation.
#' @param hyperparams Optional list of length D, each element a list with
#'                    sigma2, ell, sigma2_n (skips optimisation for that dim).
#' @return List with:
#'   \item{U_star}{n x D smoothed state (posterior mean)}
#'   \item{dU_dt}{n x D posterior mean of df/dt}
#'   \item{U_var}{n x D marginal posterior variance of f}
#'   \item{dU_var}{n x D marginal posterior variance of df/dt}
#'   \item{fits}{length-D list of gp_fit objects}
#' @export
gp_smooth <- function(U, tt, sigma2_n = NULL, hyperparams = NULL) {
  tt  <- as.vector(tt)
  n   <- nrow(U)
  D   <- ncol(U)

  sigma2_n_vec <- if (!is.null(sigma2_n)) rep_len(sigma2_n, D) else rep(NA_real_, D)

  U_star  <- matrix(0, n, D)
  dU_dt   <- matrix(0, n, D)
  U_var   <- matrix(0, n, D)
  dU_var  <- matrix(0, n, D)
  fits    <- vector("list", D)

  for (d in seq_len(D)) {
    y        <- U[, d]
    hp       <- if (!is.null(hyperparams)) hyperparams[[d]] else NULL
    sn       <- if (!is.na(sigma2_n_vec[d])) sigma2_n_vec[d] else hp$sigma2_n

    fit <- gp_fit_1d(
      tt, y,
      sigma2   = hp$sigma2,
      ell      = hp$ell,
      sigma2_n = sn
    )

    pred  <- gp_predict(fit, tt)
    dpred <- gp_predict_deriv(fit, tt)

    U_star[, d] <- pred$mean
    dU_dt[, d]  <- dpred$mean
    U_var[, d]  <- pred$var
    dU_var[, d] <- dpred$var
    fits[[d]]   <- fit
  }

  list(
    U_star = U_star,
    dU_dt  = dU_dt,
    U_var  = U_var,
    dU_var = dU_var,
    fits   = fits
  )
}
