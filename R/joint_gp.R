# WENDyGP public entry point and observed preparation. Both formulations use
# joint_gp_objective.R; latent staging lives in joint_gp_latent.R.
# Physical vectorizations are component-major (R column order).

#' Controls for WENDyGP
#'
#' @param ... Named overrides. See Details for the available controls.
#' @details Defaults are returned by calling this function without arguments.
#' Radii in gp_radius_bounds and weak_radii are in physical time units; NULL
#' selects data-based defaults. gp_nonstationary_terms is the number of smooth
#' cosine terms in the bounded log-radius function (0 gives a constant radius).
#' gp_radius_penalty penalizes their squared, frequency-weighted coefficients.
#' gp_restarts and gp_maxit control marginal-likelihood fitting. gp_quad_per_radius,
#' gp_quad_tol and gp_quad_max_refine control positive-weight convolution
#' quadrature. gp_jitter is relative numerical state-covariance stabilization.
#' kernel defaults to "matern52", using the nonstationary Matern 5/2 covariance.
#' Its radius-coefficient gradients are analytic, including the nonstationary
#' normalization factor; no finite-difference covariance evaluations are used.
#' Set gp_nonstationary_terms=0 for stationary Matern 5/2, or kernel="bump"
#' for normalized bump-convolution covariance. Weak tests remain bump functions
#' with radius choices independent of the GP kernel.
#' grid_min (default 64) selects the uniform working grid once, preserving regular observation
#' times; grid_max only caps its initial size. There is no automatic working-grid
#' refinement or optimization restart. grid_tol is the maximum
#' observation-operator interpolation RMS/noise ratio. Only observation-operator
#' accuracy can reject the grid under grid_action ("error" or "warn").
#' weak_grid_tol is a diagnostic whitened quadrature threshold, also used to
#' screen candidate weak modes; weak discrepancies do not reject a fitted grid.
#' weak_integration="grid" preserves the original grid quadrature. Experimental
#' "gp_gauss" integrates the continuous conditional GP representation, its weak
#' Jacobian and full test Gram using weak_quad_order Gauss points per fixed grid
#' interval (default 8), checked at twice that order. These are integration nodes,
#' not new optimized state variables. Exact endpoint terms are retained; em_order
#' is not applied on this path. It requires include_gp_prior=TRUE and either
#' weak_radius_method="svd" or explicit weak_radii. SVD uses quadrature-weighted
#' test values; the selected space is checked, not automatically refined/reselected.
#' weak_quad_gram_tol bounds the normalized Gram discrepancy (default 1e-3).
#' weak_design_quad_tol also checks the joint Jacobian in fixed, scaled physical
#' coordinates, independently of the optimizer's coordinate transformation.
#' Numerical failures are reported and warned about; only observation error gates
#' grid acceptance. gp_delta propagates the frozen GRID-state covariance through
#' this representation; it does not add conditional between-grid GP innovations.
#' Lagrange stencil widths are fixed within each integration interval so changes
#' in quadrature order integrate the same piecewise-polynomial representation.
#' weak_extension_tol compares GP conditional SD with fitted-state RMS. This is
#' advisory GP uncertainty, not a deterministic interpolation or quadrature error;
#' it does not trigger an accuracy warning. Final quadrature checks determine
#' quadrature_passed and weak_grid_passed; initializer checks remain in diagnostics.
#' weak_radius_method defaults to "svd": an independent dense geometric
#' multiscale pool, screened for quadrature accuracy and compressed by SVD.
#' Alternatives include "gp" (the legacy GP-radius-transfer baseline) and
#' "posterior" (an experimental posterior/noise-floor selector). The latter
#' never transfers GP radii: it evaluates a Fourier integration-error proxy using
#' its posterior mean-square plus full trajectory-covariance contribution, and
#' compares this with measurement noise propagated on the observation grid.
#' weak_radius_grid is the minimum size of its fixed diagnostic grid;
#' weak_radius_ratio spaces geometric candidates; weak_radius_tolerance is the
#' maximum posterior-proxy RMS / measurement-noise RMS (default 0.1; 1 is the
#' noise crossover). This is a fixed numerical error budget, not a calibrated
#' probability or a radius chosen using true ODE parameters.
#' This proxy does not replace the working-grid or final-solution accuracy gates.
#' Methods "svd" and experimental "sensitivity" use the same independent
#' geometric multiscale pool and screen its quadrature at a common pilot.
#' weak_design_centers=NULL places a candidate at every admissible working-grid
#' center per radius. A positive integer opts into the legacy capped placement.
#' weak_design_radii supplies physical radii or NULL uses span times
#' c(1/32,1/16,1/8,1/4,0.4). weak_design_quad_tol bounds relative quadrature
#' discrepancies in residuals and parameter sensitivities (default 1e-4).
#' weak_design_info defaults to 0.95: retain the smallest leading screened SVD
#' set reaching this fraction of the singular-value sum of the row-normalized,
#' quadrature-screened interior pool (MSG convention, not squared values).
#' The denominator precedes numerical-rank and mode-accuracy screening; an
#' unattainable target is reported, not renormalized away. weak_design_budget=NULL
#' uses this adaptive size; an integer explicitly overrides it for comparisons.
#' These are spectral information fractions, not ODE-parameter information.
#' Method "svd" retains the leading screened modes;
#' "sensitivity" reserves weak_design_coverage leading modes (default 2), then
#' greedily ranks additions by local state-adjusted parameter precision.
#' weak_design_scale supplies parameter scales, defaulting to pmax(abs(p0),1).
#' These local precision scores do not guarantee unbiased or calibrated inference.
#' weak_count caps interior centers per radius for "gp", "posterior", and
#' explicit weak_radii; it does not control the default SVD pool. bl_count is the
#' number per boundary per radius, and weak_factors is used only by method "gp".
#' Explicit weak_radii bypass automatic selection. bl_radii optionally overrides
#' boundary radii; automatic posterior selection uses only the largest pool radius
#' for the boundary block, whose EM accuracy is checked separately.
#' basis_tol and covariance_tol are relative rank tolerances.
#' ode_weighting is "test_gram" (default: quadrature-weighted L2 Gram matrix of
#' the retained tests), "gp_delta" (propagate GP trajectory covariance),
#' or "identity" (unweighted retained residual rows).
#' The Gram matrix includes interior-boundary cross terms and one identical
#' block per component; neither alternative propagates GP covariance into the
#' ODE metric. ode_units="rms_span" defaults to a dimensionless test-Gram
#' penalty: component d's precision is multiplied by T/a_d^2, where T is the
#' observed time span and a_d=sqrt(mean(Y[,d]^2)), without centering. This changes
#' the ODE penalty only, not the data or GP terms, lambda, initial pilot,
#' or SVD pool. Whitening modes are retained before this column rescaling.
#' ode_component_scale optionally replaces the RMS by a positive scalar or one
#' positive value per system component, including the latent component when
#' used with formulation="latent". These fixed scales are reported with the fit.
#' ode_units="raw" recovers the physical-unit Gram penalty. GP-delta and identity
#' always use raw units; selecting either automatically selects raw unless an
#' incompatible ode_units override is supplied. Identity depends on test-row scaling;
#' test_gram removes that dependence on the retained nonsingular subspace.
#' Weak diagnostics use the selected metric; their thresholds and lambda are
#' not numerically interchangeable across modes. GP predictive uncertainty remains
#' available for interpolation checks and state-adjusted information scoring.
#' include_bl controls the boundary block; em_order is 0, 2 or 4.
#' include_gp_prior=FALSE removes the joint GP-prior penalty (default TRUE).
#' It uses physical grid-state coordinates, with no hidden GP regularizer.
#' GP fitting still supplies initialization, optional noise estimation, and the
#' frozen gp_delta metric. The data likelihood and lambda are unchanged.
#' Without the prior the data and finite weak space may leave state directions
#' unidentified; no ridge or proper-posterior claim is substituted.
#' GP penalties have unit weight; lambda defaults to 100 and weights only the ODE term.
#' init_gls controls frozen-state GLS iterations. maxit, ftol, xtol, gtol and
#' damping control LM. maxit and gtol also set the scalar iteration budget and
#' physical stationarity tolerance; scalar native stopping tolerances remain
#' optimizer specific. bump_eta fixes the bump shape. verbose prints progress.
#' @return A named list of controls.
#' @export
wendygp_control <- function(...) {
  defaults <- list(
    kernel = "matern52",
    bump_eta = 9,
    gp_nonstationary_terms = 2L,
    gp_radius_bounds = NULL,
    gp_radius_penalty = 1,
    gp_restarts = 1L,
    gp_maxit = 200L,
    gp_quad_per_radius = 16L,
    gp_quad_tol = 1e-3,
    gp_quad_max_refine = 1L,
    gp_jitter = 1e-10,
    grid_min = 64L,
    grid_max = 500L,
    grid_tol = 0.1,
    weak_grid_tol = 0.1,
    weak_integration = "gp_gauss",
    weak_quad_order = 8L,
    weak_quad_gram_tol = 1e-3,
    weak_extension_tol = 5e-3,
    weak_extension = "lagrange",
    grid_action = "error",
    weak_radii = NULL,
    weak_radius_method = "svd",
    weak_radius_grid = 129L,
    weak_radius_ratio = sqrt(2),
    weak_radius_tolerance = 0.1,
    weak_factors = c(0.5, 1, 2),
    weak_count = 18L,
    bl_count = 1L,
    bl_radii = NULL,
    weak_design_budget = NULL,
    weak_design_coverage = 2L,
    weak_design_centers = NULL,
    weak_design_info = 0.95,
    weak_design_quad_tol = 1e-4,
    weak_design_radii = NULL,
    weak_design_scale = NULL,
    basis_tol = 1e-7,
    covariance_tol = 1e-9,
    include_bl = TRUE,
    ode_weighting = "test_gram",
    ode_units = "rms_span",
    ode_component_scale = NULL,
    include_gp_prior = TRUE,
    em_order = 4L,
    lambda = 100,
    init_gls = 5L,
    maxit = 500L,
    ftol = 1e-9,
    xtol = 1e-8,
    gtol = 1e-6,
    damping = 1e-5,
    verbose = FALSE
  )
  args <- list(...)
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))) ||
                      anyDuplicated(names(args)))) stop("Controls must have unique names.")
  if ("beta" %in% names(args))
    stop("beta has been removed. GP penalties have unit weight; use lambda for the ODE penalty.", call. = FALSE)
  unknown <- setdiff(names(args), names(defaults))
  if (length(unknown)) stop("Unknown WENDyGP control(s): ", paste(unknown, collapse = ", "))
  for (nm in names(args)) defaults[nm] <- args[nm]
  z <- defaults
  positive <- c("bump_eta", "gp_restarts", "gp_maxit", "gp_quad_per_radius",
                "gp_quad_tol", "gp_jitter", "grid_min", "grid_max", "grid_tol",
                "weak_grid_tol", "weak_count", "bl_count", "basis_tol",
                "weak_radius_grid", "weak_radius_ratio", "weak_radius_tolerance",
                "weak_design_quad_tol", "weak_quad_order", "weak_quad_gram_tol",
                "weak_extension_tol",
                "covariance_tol", "maxit", "ftol", "xtol", "gtol", "damping")
  for (nm in positive) if (!is.numeric(z[[nm]]) || length(z[[nm]]) != 1L ||
                           !is.finite(z[[nm]]) || z[[nm]] <= 0)
    stop(nm, " must be one finite positive number.")
  for (nm in c("gp_nonstationary_terms", "gp_quad_max_refine", "gp_radius_penalty",
               "init_gls", "lambda", "weak_design_coverage"))
    if (!is.numeric(z[[nm]]) || length(z[[nm]]) != 1L ||
        !is.finite(z[[nm]]) || z[[nm]] < 0) stop(nm, " must be nonnegative.")
  ints <- c("gp_nonstationary_terms", "gp_quad_max_refine", "gp_restarts", "gp_maxit",
            "grid_min", "grid_max", "weak_count", "bl_count", "init_gls", "maxit",
            "weak_radius_grid", "weak_design_coverage", "weak_quad_order")
  for (nm in ints) if (z[[nm]] != as.integer(z[[nm]])) stop(nm, " must be an integer.")
  for (nm in c("weak_design_budget", "weak_design_centers"))
    if (!is.null(z[[nm]]) && (!is.numeric(z[[nm]]) || length(z[[nm]]) != 1L ||
        !is.finite(z[[nm]]) || z[[nm]] <= 0 || z[[nm]] > .Machine$integer.max ||
        z[[nm]] != as.integer(z[[nm]]))) stop(nm, " must be NULL or a positive integer.")
  if (!is.numeric(z$weak_design_info) || length(z$weak_design_info) != 1L ||
      !is.finite(z$weak_design_info) || z$weak_design_info <= 0 || z$weak_design_info > 1)
    stop("weak_design_info must be in (0, 1].")
  if (z$grid_min < 9 || z$grid_max < z$grid_min) stop("Require 9 <= grid_min <= grid_max.")
  if (z$weak_radius_grid < 9L) stop("weak_radius_grid must be at least 9.")
  if (z$weak_radius_ratio <= 1) stop("weak_radius_ratio must exceed 1.")
  if (!is.character(z$weak_radius_method) || length(z$weak_radius_method) != 1L ||
      is.na(z$weak_radius_method) || !z$weak_radius_method %in% c("gp", "posterior", "sensitivity", "svd"))
    stop("weak_radius_method must be gp, posterior, sensitivity or svd.")
  if (!is.null(z$weak_design_budget) && z$weak_design_coverage > z$weak_design_budget)
    stop("weak_design_coverage cannot exceed weak_design_budget.")
  if (!is.character(z$ode_weighting) || length(z$ode_weighting) != 1L ||
      is.na(z$ode_weighting) || !z$ode_weighting %in% c("gp_delta", "test_gram", "identity"))
    stop("ode_weighting must be gp_delta, test_gram or identity.")
  if (!is.character(z$ode_units) || length(z$ode_units) != 1L || is.na(z$ode_units) ||
      !z$ode_units %in% c("raw", "rms_span")) stop("ode_units must be raw or rms_span.")
  if (z$ode_weighting != "test_gram") {
    if ("ode_units" %in% names(args) && z$ode_units != "raw")
      stop("ode_units='rms_span' requires ode_weighting='test_gram'.")
    z$ode_units <- "raw"
  }
  if (!z$kernel %in% c("bump", "matern52")) stop("kernel must be bump or matern52.")
  if (!z$grid_action %in% c("error", "warn")) stop("grid_action must be error or warn.")
  if (length(z$em_order) != 1L || !z$em_order %in% c(0L, 2L, 4L)) stop("em_order must be 0, 2 or 4.")
  for (nm in c("include_bl", "include_gp_prior", "verbose"))
    if (!is.logical(z[[nm]]) || length(z[[nm]]) != 1L || is.na(z[[nm]])) stop(nm, " must be logical.")
  if (!is.character(z$weak_extension) || !length(z$weak_extension) ||
      anyNA(z$weak_extension) || !all(z$weak_extension %in% c("lagrange", "gp")))
    stop("weak_extension must be lagrange or gp, one entry or one per component.")
  if (!is.character(z$weak_integration) || length(z$weak_integration) != 1L ||
      is.na(z$weak_integration) || !z$weak_integration %in% c("grid", "gp_gauss"))
    stop("weak_integration must be grid or gp_gauss.")
  if (z$weak_quad_order < 2L || z$weak_quad_order > 32L)
    stop("weak_quad_order must be an integer from 2 to 32.")
  # gp_gauss serves neither a free state nor a non-SVD radius method, so such a
  # request is only answerable by grid quadrature. Controls here are routinely
  # built complete and then mutated field by field, and solveWendyGP re-resolves
  # them, so whether the caller named weak_integration is not recoverable at
  # this point -- a precedence rule on names(args) cannot work. Coerce to the
  # only scheme that can answer, and say so rather than doing it silently.
  if (z$weak_integration == "gp_gauss" && (!z$include_gp_prior ||
      (is.null(z$weak_radii) && z$weak_radius_method != "svd"))) {
    z$weak_integration <- "grid"
    warning("weak_integration coerced to grid: gp_gauss supports neither a free ",
            "state (include_gp_prior=FALSE) nor weak_radius_method other than ",
            "svd without explicit weak_radii.", call. = FALSE)
  }
  for (nm in c("weak_radii", "weak_factors", "gp_radius_bounds", "bl_radii",
               "weak_design_radii", "weak_design_scale", "ode_component_scale"))
    if (!is.null(z[[nm]]) && (!is.numeric(z[[nm]]) || !length(z[[nm]]) ||
                            any(!is.finite(z[[nm]]) | z[[nm]] <= 0))) stop(nm, " must be positive.")
  if (!is.null(z$gp_radius_bounds) && (length(z$gp_radius_bounds) != 2L ||
                                    diff(z$gp_radius_bounds) <= 0)) stop("Invalid gp_radius_bounds.")
  z
}

# The peak rescaling exp(eta) cancels under L2 normalization but avoids tiny
# numbers. The derivatives come from the SAME phi and derivative cache used
# by the existing weak-form code. Evaluation outside support is exactly zero.
.jgp_bump <- function(x, order = 0L, eta = 9) {
  out <- x * 0
  inside <- abs(x) < 1 & is.finite(x)
  if (any(inside)) {
    fun <- function(t, r) phi(t, r, eta = eta)
    deriv <- test_function_derivative(fun, 1, 1, order)
    out[inside] <- exp(eta) * deriv(x[inside])
  }
  out
}

.jgp_chol <- function(K, jitter = 0, label = "covariance") {
  K <- (K + t(K)) / 2
  added <- jitter * max(diag(K))
  R <- tryCatch(chol(K + diag(added, nrow(K))), error = function(e) NULL)
  if (is.null(R)) stop(label, " is not positive definite at the requested stabilization.")
  list(R = R, added = added)
}

.jgp_solve <- function(R, b) backsolve(R, forwardsolve(t(R), b))

.jgp_radius <- function(x, a, bounds) {
  # Preserve the location-by-coefficient Jacobian shape even for one location.
  Z <- matrix(vapply(0:(length(a) - 1L), function(j) cos(j * pi * x),
                     numeric(length(x))), nrow = length(x), ncol = length(a))
  s <- stats::plogis(as.vector(Z %*% a))
  width <- diff(log(bounds))
  list(r = exp(log(bounds[1]) + width * s),
       derivative = Z * as.vector(width * s * (1 - s)))
}

.jgp_quad <- function(bounds, per_radius) {
  # Integration domain contains every admissible support for x in [0,1].
  lo <- -bounds[2]; hi <- 1 + bounds[2]
  intervals <- ceiling((hi - lo) / (bounds[1] / per_radius))
  nodes <- seq(lo, hi, length.out = intervals + 1L)
  w <- rep((hi - lo) / intervals, length(nodes))
  w[c(1L, length(w))] <- w[c(1L, length(w))] / 2
  list(nodes = nodes, weights = w)
}

.jgp_features <- function(x, a, bounds, quad, eta, derivatives = FALSE) {
  rad <- .jgp_radius(x, a, bounds)
  v <- sweep(outer(x, quad$nodes, "-"), 1, rad$r, "/")
  b <- .jgp_bump(v, eta = eta)
  raw <- sweep(b, 2, sqrt(quad$weights), "*") / sqrt(rad$r)
  norms <- sqrt(rowSums(raw^2))
  if (any(!is.finite(norms) | norms <= 0)) stop("Unresolved convolution quadrature.")
  B <- raw / norms
  dB <- NULL
  if (derivatives) {
    raw_deriv <- sweep(-v * .jgp_bump(v, 1L, eta) - 0.5 * b,
                      2, sqrt(quad$weights), "*") / sqrt(rad$r)
    dB <- (raw_deriv - B * rowSums(B * raw_deriv)) / norms
  }
  list(B = B, dB = dB, radius = rad$r, dr = rad$derivative)
}

.jgp_correlation <- function(x, a, bounds, quad, control, derivatives = FALSE) {
  if (control$kernel == "bump") {
    feat <- .jgp_features(x, a, bounds, quad, control$bump_eta, derivatives)
    R <- tcrossprod(feat$B)
    dR <- if (derivatives) lapply(seq_along(a), function(j) {
      part <- tcrossprod(feat$dB * feat$dr[, j], feat$B)
      part + t(part)
    }) else NULL
    return(list(R = R, dR = dR, features = feat$B))
  }
  # Paciorek-Schervish determinant factor and averaged squared local scales.
  radius <- .jgp_radius(x, a, bounds)
  rad <- radius$r; rad2 <- rad^2
  avg <- outer(rad2, rad2, "+") / 2
  z <- sqrt(5) * abs(outer(x, x, "-")) / sqrt(avg)
  prefactor <- sqrt(outer(rad, rad) / avg); decay <- exp(-z)
  R <- prefactor * (1 + z + z^2 / 3) * decay
  # radius$derivative is d log(r_i) / d a_j. Differentiate both the
  # determinant prefactor and the Matern factor through the averaged r^2.
  # This form never divides by distance or z, so repeated times and diagonal
  # entries are well defined (the unit-diagonal correlation has zero gradient).
  radial <- if (derivatives) prefactor * decay * z^2 * (1 + z) / 6
  dR <- if (derivatives) lapply(seq_along(a), function(j) {
    b <- radius$derivative[, j]
    dlogavg <- outer(rad2 * b, rad2 * b, "+") / avg
    R * (outer(b, b, "+") - dlogavg) / 2 + radial * dlogavg
  }) else NULL
  list(R = R, dR = dR, features = NULL)
}

.jgp_gp_fit <- function(tt, y, noise_sd, control) {
  origin <- min(tt); span <- diff(range(tt)); x <- (tt - origin) / span
  n <- length(y); center <- mean(y); scale <- stats::sd(y)
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps) * max(1, abs(center)))
    stop("A component has essentially constant observations; its GP scale is not identifiable.")
  ys <- (y - center) / scale
  bounds <- if (is.null(control$gp_radius_bounds))
    c(max(0.02, min(0.2, 2 * stats::median(diff(x)))), 2) else control$gp_radius_bounds / span
  q <- control$gp_nonstationary_terms + 1L
  penalty_weights <- (0:(q - 1L))^4 * control$gp_radius_penalty
  known <- !is.null(noise_sd)
  noise2 <- if (known) (noise_sd / scale)^2 else NULL
  quad <- .jgp_quad(bounds, control$gp_quad_per_radius)
  best <- NULL
  for (refine in 0:control$gp_quad_max_refine) {
    cached_par <- cached <- NULL
    evaluate <- function(par) {
      if (identical(par, cached_par)) return(cached)
      a <- par[seq_len(q)]; s <- exp(par[q + 1L])
      cr <- .jgp_correlation(x, a, bounds, quad, control, TRUE)
      A <- if (known) s * cr$R + diag(noise2, n) else cr$R + diag(s, n)
      R <- tryCatch(chol(A), error = function(e) NULL)
      if (is.null(R)) return(list(value = 1e50, gradient = rep(0, length(par))))
      Ai <- .jgp_solve(R, diag(n)); ai1 <- rowSums(Ai)
      mu <- sum(ai1 * ys) / sum(ai1); e <- ys - mu
      alpha <- as.vector(Ai %*% e)
      tau2 <- if (known) s else max(sum(e * alpha) / n, .Machine$double.eps)
      value <- sum(log(diag(R))) + if (known) sum(e * alpha) / 2 else n * log(tau2) / 2
      value <- value + sum(penalty_weights * a^2) / 2
      score <- (Ai - tcrossprod(alpha) / if (known) 1 else tau2) / 2
      gradient <- vapply(cr$dR, function(dR) sum(score * dR) * if (known) s else 1, numeric(1))
      gradient <- c(gradient + penalty_weights * a,
                    if (known) sum(score * (s * cr$R)) else s * sum(diag(score)))
      cached_par <<- par
      cached <<- list(value = value, gradient = gradient, mean = mu, tau2 = tau2,
                     noise2 = if (known) noise2 else tau2 * s, R = cr$R)
      cached
    }
    starts <- if (!is.null(best)) list(best$par) else lapply(seq_len(control$gp_restarts), function(j) {
      frac <- seq(0.3, 0.8, length.out = control$gp_restarts)[j]
      c(stats::qlogis(frac), rep(0, q - 1L), if (known) 0 else log(0.03))
    })
    fits <- lapply(starts, function(init) stats::optim(
      init, function(a) evaluate(a)$value, function(a) evaluate(a)$gradient,
      method = "L-BFGS-B", lower = c(rep(-8, q), -18), upper = c(rep(8, q), 10),
      control = list(maxit = control$gp_maxit, factr = 1e7)))
    best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
    val <- evaluate(best$par)
    finer <- .jgp_quad(bounds, control$gp_quad_per_radius * 2^(refine + 1L))
    fine_R <- .jgp_correlation(x, best$par[seq_len(q)], bounds, finer, control)$R
    quad_error <- max(abs(val$R - fine_R))
    if (control$kernel != "bump" || quad_error <= control$gp_quad_tol) break
    if (refine == control$gp_quad_max_refine) stop("GP convolution quadrature did not converge.")
    quad <- finer
  }
  Kyy <- val$tau2 * scale^2 * val$R + diag(val$noise2 * scale^2, n)
  Ryy <- .jgp_chol(Kyy, label = "GP observation covariance")$R
  structure(list(tt = tt, x = x, y = y, origin = origin, span = span,
    mean = center + scale * val$mean, tau2 = scale^2 * val$tau2,
    noise2 = scale^2 * val$noise2, radius_coef = best$par[seq_len(q)],
    radius_bounds = bounds, quad = quad, control = control, Ryy = Ryy,
    alpha = .jgp_solve(Ryy, y - center - scale * val$mean),
    convergence = best$convergence, message = best$message, objective = best$value,
    quadrature_error = quad_error, quadrature_refinements = refine), class = "wendygp_gp")
}

.jgp_gp_predict <- function(fit, tt) {
  # Evaluate training and prediction locations with the same feature quadrature.
  x <- (tt - fit$origin) / fit$span
  allx <- c(fit$x, x); n <- length(fit$x); ix <- n + seq_along(x)
  R <- .jgp_correlation(allx, fit$radius_coef, fit$radius_bounds, fit$quad, fit$control)$R
  Ks <- fit$tau2 * R[ix, seq_len(n), drop = FALSE]
  K <- fit$tau2 * R[ix, ix, drop = FALSE]
  V <- forwardsolve(t(fit$Ryy), t(Ks))
  Sigma <- K - crossprod(V)
  list(mean = as.vector(fit$mean + Ks %*% fit$alpha), K = K,
       Sigma = (Sigma + t(Sigma)) / 2,
       radius = fit$span * .jgp_radius(x, fit$radius_coef, fit$radius_bounds)$r)
}

.jgp_gp_cross <- function(fit, t1, t2) {
  x <- (t1 - fit$origin) / fit$span; y <- (t2 - fit$origin) / fit$span
  if (fit$control$kernel == "bump") {
    Bx <- .jgp_features(x, fit$radius_coef, fit$radius_bounds, fit$quad, fit$control$bump_eta)$B
    By <- .jgp_features(y, fit$radius_coef, fit$radius_bounds, fit$quad, fit$control$bump_eta)$B
    return(fit$tau2 * tcrossprod(Bx, By))
  }
  rx <- .jgp_radius(x, fit$radius_coef, fit$radius_bounds)$r
  ry <- .jgp_radius(y, fit$radius_coef, fit$radius_bounds)$r
  avg <- outer(rx^2, ry^2, "+") / 2
  z <- sqrt(5) * abs(outer(x, y, "-")) / sqrt(avg)
  fit$tau2 * sqrt(outer(rx, ry) / avg) * (1 + z + z^2 / 3) * exp(-z)
}

.jgp_model_cache <- new.env(parent = emptyenv())

.jgp_model <- function(f, D, J, em_order) {
  u <- sym_build(lapply(seq_len(D), function(j) sym_symbol(paste0("u", j))))
  p <- sym_build(lapply(seq_len(J), function(j) sym_symbol(paste0("p", j))))
  time <- sym_symbol("t"); vars <- c(p, u, time)
  rhs <- f(u, p, time)
  if (sym_length(rhs) != D) stop("f must return one derivative per observed component.")
  # Compare evaluated expressions, not only closure identity: captured numeric
  # constants may have changed in the caller's environment since the last fit.
  signature <- sym_strings(rhs)
  old <- .jgp_model_cache$last
  if (!is.null(old) && identical(old$signature, signature) && old$D == D && old$J == J &&
      old$em_order == em_order && identical(old$backend, sym_backend())) return(old)
  jets <- list(rhs)
  njet <- if (em_order == 4L) 4L else if (em_order == 2L) 2L else 1L
  if (njet > 1L) for (j in 2:njet)
    jets[[j]] <- compute_symbolic_total_time_deriv(jets[[j - 1L]], u, rhs, time)
  px <- c(p, u)
  # Sparsity and constancy of df/du, decided symbolically once per model. A
  # structurally zero block contributes nothing to the state Jacobian; a block
  # free of both u and t is the SAME scalar at every quadrature node, so its
  # weak product collapses to a scalar multiple of a matrix that is fixed for
  # the whole solve. Row order of sym_strings is column-major over the D-by-D
  # Jacobian, so both masks index as [component, state].
  ju_sym <- compute_symbolic_jacobian(rhs, u)
  ju_str <- sym_strings(ju_sym)
  ju_dep <- matrix(sym_strings(compute_symbolic_jacobian(ju_sym, c(u, time))), D * D, D + 1L)
  model <- list(f_user = f, signature = signature, D = D, J = J, em_order = em_order, backend = sym_backend(),
    jet = lapply(jets, build_fn, vars = vars),
    jet_jac = lapply(jets, function(z) build_fn(compute_symbolic_jacobian(z, px), vars)),
    ju_zero = matrix(ju_str == "0", D, D),
    ju_const = matrix(apply(ju_dep == "0", 1L, all), D, D),
    affine = all(sym_strings(compute_symbolic_jacobian(
      compute_symbolic_jacobian(rhs, p), p)) == "0"))
  .jgp_model_cache$last <- model
  model
}

# Physical trapezoidal weights, including irregular observation grids. The
# diagnostic grid and observation grid must use the SAME integral convention.
.jgp_time_weights <- function(tt) {
  h <- diff(tt)
  c(h[1], head(h, -1L) + tail(h, -1L), h[length(h)]) / 2
}

# E ||A U||^2 / nrow(A), conditional on fitted GP hyperparameters. B maps
# ORIGINAL measurement errors, not independent pseudo-observations on a dense
# grid. No numerical GP jitter is added to either stochastic quantity here.
.jgp_proxy_moments <- function(A, mean, Sigma, B, noise2) {
  mean_square <- mean(as.vector(A %*% mean)^2)
  terms <- (A %*% Sigma) * A
  posterior_variance <- sum(terms) / nrow(A)
  roundoff <- 1e-8 * max(abs(diag(Sigma))) * sum(A^2) / nrow(A)
  if (posterior_variance < -max(roundoff, .Machine$double.eps))
    stop("Posterior radius diagnostic encountered a negative covariance contribution.")
  posterior_variance <- max(0, posterior_variance)
  noise_variance <- noise2 * sum(B^2) / nrow(B)
  total <- mean_square + posterior_variance
  list(mean_square = mean_square, posterior_variance = posterior_variance,
       expected_square = total, noise_variance = noise_variance,
       rms_ratio = if (noise_variance > 0) sqrt(total / noise_variance) else Inf)
}

.jgp_radius_centers <- function(tt, radius, control, dense = FALSE) {
  span <- diff(range(tt))
  count <- max(3L, min(control$weak_count,
    ceiling((span - 2 * radius) / (radius / if (dense) 4 else 2)) + 1L))
  seq(min(tt) + radius, max(tt) - radius, length.out = count)
}

# Posterior version of the high-frequency proxy used by
# find_min_radius_int_error. Both Fourier phases are retained to avoid a zero
# caused merely by phase alignment. Its carrier frequency is tied to the
# ORIGINAL observation count, not the number of interpolated grid values.
.jgp_radius_operators <- function(tt, diagnostic_tt, radius, control) {
  centers <- .jgp_radius_centers(tt, radius, control, dense = TRUE)
  design <- data.frame(center = centers, radius = radius)
  w <- .jgp_time_weights(diagnostic_tt); wo <- .jgp_time_weights(tt)
  raw <- .jgp_test_rows(diagnostic_tt, design, 0L, control$bump_eta)
  norms <- sqrt(as.vector(raw^2 %*% w))
  if (any(!is.finite(norms) | norms <= 0)) stop("Unresolved radius diagnostic supports.")
  wave <- max(1L, floor(length(tt) / 3L) - 1L)
  phase <- function(t) 2 * pi * wave * (t - min(tt)) / diff(range(tt))
  operator <- function(values, times, weights) {
    values <- values / norms
    rbind(sweep(values, 2, weights * sin(phase(times)), "*"),
          sweep(values, 2, weights * cos(phase(times)), "*"))
  }
  list(A = operator(raw, diagnostic_tt, w),
       B = operator(.jgp_test_rows(tt, design, 0L, control$bump_eta), tt, wo),
       centers = centers, wavenumber = wave)
}

.jgp_posterior_radius <- function(tt, fits, control) {
  # A dedicated diagnostic grid makes physical test supports invariant to
  # working-grid choice. Its resolution can be checked independently.
  size <- max(length(tt), control$weak_radius_grid)
  diagnostic_tt <- seq(min(tt), max(tt), length.out = size)
  span <- diff(range(tt)); smallest <- 4 * span / (size - 1L); largest <- 0.4 * span
  smallest <- min(smallest, largest)
  steps <- floor(log(largest / smallest) / log(control$weak_radius_ratio))
  radii <- sort(unique(c(smallest * control$weak_radius_ratio^(0:steps), largest)))
  predictions <- lapply(fits, .jgp_gp_predict, tt = diagnostic_tt)
  operators <- lapply(radii, function(r) .jgp_radius_operators(tt, diagnostic_tt, r, control))
  records <- list(); index <- 0L
  for (d in seq_along(fits)) for (j in seq_along(radii)) {
    p <- predictions[[d]]; op <- operators[[j]]
    moments <- .jgp_proxy_moments(op$A, p$mean, p$Sigma, op$B, fits[[d]]$noise2)
    index <- index + 1L
    records[[index]] <- data.frame(component = d, radius = radii[j],
      mean_square = moments$mean_square, posterior_variance = moments$posterior_variance,
      expected_square = moments$expected_square, noise_variance = moments$noise_variance,
      rms_ratio = moments$rms_ratio, centers = length(op$centers),
      passed = is.finite(moments$rms_ratio) && moments$rms_ratio <= control$weak_radius_tolerance)
  }
  profile <- do.call(rbind, records)
  per_component <- vapply(seq_along(fits), function(d) {
    rows <- profile[profile$component == d, ]
    # An isolated Fourier notch must not admit an otherwise unresolved family.
    ok <- which(as.logical(rev(cumprod(rev(rows$passed)))))
    if (length(ok)) rows$radius[ok[1L]] else NA_real_
  }, numeric(1))
  selection <- list(method = "posterior", profile = profile, candidates = radii,
    component_minimum = per_component, diagnostic_tt = diagnostic_tt,
    observation_points = length(tt), wavenumber = operators[[1L]]$wavenumber,
    tolerance = control$weak_radius_tolerance, includes_posterior_covariance = TRUE,
    noise_on_observation_grid = TRUE, hyperparameters_conditioned_on = TRUE,
    criterion = "posterior Fourier-proxy RMS / original-observation noise RMS",
    passed = all(is.finite(per_component)))
  if (!selection$passed) {
    bad <- which(!is.finite(per_component))
    condition <- structure(list(message = paste0("No posterior/noise-floor weak radius is admissible for component(s) ",
      paste(bad, collapse = ", "), ". Inspect the radius profile, GP/noise fit and observation resolution; ",
      "no radius was silently substituted."), call = NULL, radius_selection = selection),
      class = c("wendygp_radius_error", "error", "condition"))
    stop(condition)
  }
  # The current residual shares one test bank across components, so require the
  # entire retained family to pass for EVERY component, not only a pooled average.
  selection$minimum <- max(per_component)
  selection$radii <- radii[radii >= selection$minimum]
  selection
}

.jgp_test_design <- function(tt, fits, control) {
  span <- diff(range(tt)); h <- min(diff(tt))
  radii <- control$weak_radii
  automatic <- is.null(radii)
  selection <- list(method = if (automatic) control$weak_radius_method else "explicit")
  if (automatic && control$weak_radius_method == "posterior") {
    selection <- .jgp_posterior_radius(tt, fits, control)
    radii <- selection$radii
  } else if (automatic) {
    gp_radii <- vapply(fits, function(g) g$span * stats::median(
      .jgp_radius(seq(0, 1, length.out = 25), g$radius_coef, g$radius_bounds)$r), numeric(1))
    radii <- as.vector(outer(gp_radii, control$weak_factors))
    # Physical supports do not change when the working grid is refined.
    radii <- pmin(0.4 * span, pmax(4 * h, radii))
    # Snap proposals to a shared geometric family. Nearly identical radii from
    # different components must not create a forest of derivative-like BL modes.
    family <- sort(unique(pmin(0.4 * span, pmax(4 * h, span * 2^(-6:-1)))))
    radii <- sort(unique(vapply(radii, function(r) family[which.min(abs(log(family / r)))], numeric(1))))
  }
  radii <- sort(unique(radii))
  if (any(radii >= span / 2)) stop("weak_radii must be smaller than half the time span.")
  a <- min(tt); b <- max(tt)
  interior <- do.call(rbind, lapply(radii, function(r) {
    data.frame(center = .jgp_radius_centers(tt, r, control,
      dense = automatic && control$weak_radius_method == "posterior"), radius = r)
  }))
  boundary_radii <- control$bl_radii
  if (is.null(boundary_radii)) boundary_radii <-
    if (automatic && control$weak_radius_method == "posterior") max(radii) else radii
  if (any(boundary_radii >= span / 2)) stop("bl_radii must be smaller than half the time span.")
  boundary_radii <- sort(unique(boundary_radii))
  boundary <- if (control$include_bl) do.call(rbind, lapply(boundary_radii, function(r) {
    offset <- seq(0, 0.4 * r, length.out = control$bl_count)
    data.frame(center = c(a + offset, b - offset), radius = r)
  })) else data.frame(center = numeric(), radius = numeric())
  list(interior = interior, boundary = boundary, radii = radii,
       radius_selection = selection, boundary_radii = boundary_radii)
}

.jgp_test_rows <- function(tt, design, order, eta) {
  if (!nrow(design)) return(matrix(0, 0, length(tt)))
  v <- (matrix(rep(tt,each=nrow(design)),nrow=nrow(design))-design$center)/design$radius
  .jgp_bump(v, order, eta) / design$radius^order
}

.jgp_basis_map <- function(V, tol, spectrum = FALSE) {
  if (!nrow(V)) return(matrix(0, 0, 0))
  norms <- sqrt(rowSums(V^2))
  if (any(norms == 0)) stop("Weak supports are unresolved on the initial grid.")
  X <- V / norms
  decomp <- svd(X, nu = min(dim(X)), nv = 0)
  keep <- which(decomp$d > max(decomp$d) * tol)
  map <- sweep(t(decomp$u[, keep, drop = FALSE]) / decomp$d[keep], 2, norms, "/")
  if (spectrum) list(map = map, singular_values = decomp$d,
                     retained = keep) else map
}

.jgp_weak <- function(tt, design, control, maps = NULL) {
  m <- length(tt); h <- mean(diff(tt))
  weights <- rep(h, m); weights[c(1, m)] <- h / 2
  ints <- .jgp_test_rows(tt, design$interior, 0L, control$bump_eta)
  bls <- .jgp_test_rows(tt, design$boundary, 0L, control$bump_eta)
  if (is.null(maps)) maps <- list(
    interior = .jgp_basis_map(ints, control$basis_tol),
    boundary = .jgp_basis_map(bls, control$basis_tol))
  V <- rbind(maps$interior %*% ints, maps$boundary %*% bls)
  Vp <- rbind(maps$interior %*% .jgp_test_rows(tt, design$interior, 1L, control$bump_eta),
              maps$boundary %*% .jgp_test_rows(tt, design$boundary, 1L, control$bump_eta))
  ni <- nrow(maps$interior); nb <- nrow(maps$boundary)
  endpoints <- lapply(c(tt[1], tt[m]), function(t) {
    raw <- vapply(0:4, function(j) as.vector(.jgp_test_rows(
      t, design$boundary, j, control$bump_eta)), numeric(nrow(design$boundary)))
    if (!nb) matrix(0, 0, 5) else maps$boundary %*% raw
  })
  V <- sweep(V, 2, weights, "*"); Vp <- sweep(Vp, 2, weights, "*")
  if (nb) {
    Vp[ni + seq_len(nb), 1] <- Vp[ni + seq_len(nb), 1] + endpoints[[1]][, 1]
    Vp[ni + seq_len(nb), m] <- Vp[ni + seq_len(nb), m] - endpoints[[2]][, 1]
  }
  list(V = V, Vp = Vp, ni = ni, nb = nb, K = ni + nb, tt = tt,
       endpoints = endpoints, maps = maps, design = design, h = h,
       em_order = control$em_order)
}

.jgp_weak_eval <- function(U, p, weak, model, jacobian = TRUE, whiten = FALSE,
                           adjoint = FALSE) {
  if (identical(weak$integration, "gp_gauss"))
    return(.jgp_gauss_eval(U, p, weak, model, jacobian, whiten, adjoint))
  if (whiten) stop("Whitened state Jacobians are only built by gp_gauss integration.")
  m <- nrow(U); D <- model$D; J <- model$J; K <- weak$K
  input <- rbind(matrix(p, J, m), t(U), weak$tt)
  F <- model$jet[[1]](input)
  r <- weak$V %*% F + weak$Vp %*% U
  Jp <- Ju <- NULL
  if (jacobian) {
    df <- array(model$jet_jac[[1]](input), c(m, D, J + D))
    Jp <- matrix(weak$V %*% matrix(df[, , seq_len(J), drop = FALSE], m, D * J), K * D, J)
    Ju <- matrix(0, K * D, m * D)
    for (a in seq_len(D)) for (b in seq_len(D)) {
      ir <- (a - 1L) * K + seq_len(K); ic <- (b - 1L) * m + seq_len(m)
      Ju[ir, ic] <- sweep(weak$V, 2, df[, a, J + b], "*") + if (a == b) weak$Vp else 0
    }
  }
  if (weak$nb && weak$em_order) {
    brows <- weak$ni + seq_len(weak$nb)
    for (side in 1:2) {
      it <- if (side == 1L) 1L else m
      sgn <- if (side == 1L) -1 else 1
      inp <- matrix(c(p, U[it, ], weak$tt[it]), ncol = 1L)
      fd <- lapply(model$jet, function(fn) as.vector(fn(inp)))
      coeff <- matrix(0, 5L, D)
      coeff[1:3, ] <- -weak$h^2 / 12 * g_coeffs(fd, U[it, ], 1L)
      if (weak$em_order == 4L) coeff <- coeff + weak$h^4 / 720 * g_coeffs(fd, U[it, ], 3L)
      r[brows, ] <- r[brows, ] + sgn * weak$endpoints[[side]] %*% coeff
      if (jacobian) {
        fdj <- lapply(model$jet_jac, function(fn) as.vector(fn(inp)))
        du <- as.vector(cbind(matrix(0, D, J), diag(D)))
        coefj <- matrix(0, 5L, D * (J + D))
        coefj[1:3, ] <- -weak$h^2 / 12 * g_coeffs(fdj, du, 1L)
        if (weak$em_order == 4L) coefj <- coefj + weak$h^4 / 720 * g_coeffs(fdj, du, 3L)
        ej <- array(sgn * weak$endpoints[[side]] %*% coefj, c(weak$nb, D, J + D))
        for (a in seq_len(D)) {
          ir <- (a - 1L) * K + brows
          Jp[ir, ] <- Jp[ir, , drop = FALSE] + matrix(ej[, a, seq_len(J), drop = FALSE], weak$nb, J)
          for (b in seq_len(D)) {
            ic <- (b - 1L) * m + it
            Ju[ir, ic] <- Ju[ir, ic] + ej[, a, J + b]
          }
        }
      }
    }
  }
  list(r = as.vector(r), Jp = Jp, Ju = Ju)
}

# Rank decisions are made on the correlation matrix. They are invariant under
# independent nonzero row rescalings, including sign flips. Discarded modes
# define an explicitly reported projected penalty, never a hidden raw ridge.
.jgp_whitener <- function(S, tol) {
  if (!nrow(S)) return(list(W = matrix(0, 0, 0), rank = 0L, discarded = 0L))
  S <- (S + t(S)) / 2
  sd <- sqrt(pmax(diag(S), 0))
  active <- which(sd > 0)
  if (!length(active)) stop("Residual covariance has zero rank.")
  corr <- S[active, active, drop = FALSE] / outer(sd[active], sd[active])
  e <- eigen(corr, symmetric = TRUE)
  if (min(e$values) < -1e-6 * max(e$values)) stop("Residual covariance is not positive semidefinite.")
  keep <- which(e$values > max(e$values) * tol)
  W <- matrix(0, length(keep), nrow(S))
  W[, active] <- sweep(t(e$vectors[, keep, drop = FALSE]) / sqrt(e$values[keep]),
                        2, sd[active], "/")
  list(W = W, rank = length(keep), discarded = nrow(S) - length(keep),
       smallest_retained = min(e$values[keep]), tolerance = tol)
}

.jgp_weights <- function(S, weak, D, tol) {
  ii <- unlist(lapply(seq_len(D), function(d) (d - 1L) * weak$K + seq_len(weak$ni)))
  bb <- setdiff(seq_len(nrow(S)), ii)
  wi <- .jgp_whitener(S[ii, ii, drop = FALSE], tol)
  Wi <- matrix(0, wi$rank, nrow(S)); Wi[, ii] <- wi$W
  if (length(bb)) {
    C <- S[ii, bb, drop = FALSE]
    B <- wi$W %*% C
    conditional <- S[bb, bb, drop = FALSE] - crossprod(B)
    # Tiny negative roundoff on a Schur diagonal is not a covariance ridge.
    diag(conditional) <- pmax(diag(conditional), 0)
    wb <- .jgp_whitener(conditional, tol)
    T <- matrix(0, length(bb), nrow(S)); T[, bb] <- diag(length(bb))
    T[, ii] <- -crossprod(B, wi$W)
    Wb <- wb$W %*% T
  } else {
    wb <- list(rank = 0L, discarded = 0L); conditional <- matrix(0, 0, 0)
    Wb <- matrix(0, 0, nrow(S))
  }
  list(W = rbind(Wi, Wb), Wi = Wi, Wb = Wb, interior = wi, boundary = wb,
       conditional_covariance = conditional, ii = ii, bb = bb)
}

.jgp_propagate <- function(Ju, Sigma, m, D) {
  S <- matrix(0, nrow(Ju), nrow(Ju))
  for (d in seq_len(D)) {
    j <- Ju[, (d - 1L) * m + seq_len(m), drop = FALSE]
    S <- S + j %*% Sigma[[d]] %*% t(j)
  }
  (S + t(S)) / 2
}

# Private helpers also accept saved pre-switch control lists.
.jgp_ode_mode <- function(control) {
  if (is.null(control$ode_weighting)) "gp_delta" else control$ode_weighting
}

# Missing fields in archived controls retain the original raw-unit objective.
.jgp_ode_scaling <- function(Y, tt, control) {
  units <- if (is.null(control$ode_units)) "raw" else control$ode_units
  if (.jgp_ode_mode(control) != "test_gram") units <- "raw"
  rms <- sqrt(colMeans(Y^2)); span <- diff(range(tt))
  if (!is.null(control$ode_component_scale)) {
    if (!length(control$ode_component_scale) %in% c(1L, ncol(Y)))
      stop("ode_component_scale must be scalar or have one value per component.")
    rms <- rep_len(control$ode_component_scale, ncol(Y))
  }
  q <- rep(1, ncol(Y))
  if (units == "rms_span") {
    if (!is.finite(span) || span <= 0 || any(!is.finite(rms) | rms <= 0))
      stop("RMS/span scaling requires a positive span and positive finite component RMS.")
    q <- span / rms^2
    if (any(!is.finite(q) | q <= 0)) stop("Nonfinite RMS/span precision; rescale input units.")
  }
  list(units = units, rms = rms, span = span, precision = q,
       scale_source = if (is.null(control$ode_component_scale)) "observed_rms" else "supplied",
       centered = FALSE, frozen = TRUE)
}

# Scale the original retained whiteners, never re-truncate a unit-scaled Gram.
.jgp_scale_weights <- function(weights, q, K) {
  stopifnot(all(is.finite(q) & q > 0), ncol(weights$W) == K * length(q))
  if (all(q == 1)) return(weights)
  factors <- rep(sqrt(q), each = K); out <- weights
  for (k in c("W", "Wi", "Wb")) out[[k]] <- sweep(weights[[k]], 2, factors, "*")
  out$interior$W <- sweep(weights$interior$W, 2, factors[weights$ii], "*")
  if (length(weights$bb)) {
    out$boundary$W <- sweep(weights$boundary$W, 2, factors[weights$bb], "*")
    out$conditional_covariance <- weights$conditional_covariance /
      outer(factors[weights$bb], factors[weights$bb])
  }
  out
}

.jgp_fit_metric <- function(metric, weak, Y, tt, control) {
  scaling <- .jgp_ode_scaling(Y, tt, control)
  weights <- .jgp_ode_weights(metric, weak, ncol(Y), control)
  factors <- rep(sqrt(scaling$precision), each = weak$K)
  list(Omega = metric / outer(factors, factors), raw = metric,
       weights = .jgp_scale_weights(weights, scaling$precision, weak$K), scaling = scaling)
}

# The matrix traditionally stored as Omega is a residual METRIC for the two
# non-GP modes, not an estimate of residual sampling covariance. All rows use
# the existing component-major order. Optional rows avoids full propagation
# during interior-only initialization and gives consistent principal blocks.
.jgp_ode_metric <- function(weak, D, control, Ju = NULL, Sigma = NULL, rows = NULL) {
  mode <- .jgp_ode_mode(control)
  if (is.null(rows)) rows <- seq_len(weak$K * D)
  if (mode == "gp_delta")
    return(.jgp_propagate(Ju[rows, , drop = FALSE], Sigma, length(weak$tt), D))
  if (mode == "identity") return(diag(length(rows)))
  if (mode != "test_gram") stop("Unknown ODE weighting mode.")
  if (identical(weak$integration, "gp_gauss")) {
    metric <- kronecker(diag(D), weak$gram)
    return(metric[rows, rows, drop = FALSE])
  }
  # V already contains trapezoidal weights: V V' would incorrectly use dt^2.
  # Recover Phi sqrt(w), giving integral Phi_a(t) Phi_b(t) dt, with boundary
  # tests included. Endpoint/EM derivative corrections belong to the residual,
  # not to this test-function geometry.
  w <- rep(weak$h, length(weak$tt)); w[c(1L, length(w))] <- weak$h / 2
  M <- tcrossprod(sweep(weak$V, 2, sqrt(w), "/"))
  metric <- kronecker(diag(D), M)
  metric[rows, rows, drop = FALSE]
}

.jgp_ode_weights <- function(metric, weak, D, control) {
  if (.jgp_ode_mode(control) != "identity")
    return(.jgp_weights(metric, weak, D, control$covariance_tol))
  # No eigenanalysis, row rescaling, or rank truncation in ordinary LS.
  ii <- unlist(lapply(seq_len(D), function(d) (d - 1L) * weak$K + seq_len(weak$ni)))
  bb <- setdiff(seq_len(weak$K * D), ii)
  I <- diag(weak$K * D)
  Wi <- I[ii, , drop = FALSE]; Wb <- I[bb, , drop = FALSE]
  info <- function(n) list(rank = n, discarded = 0L,
    smallest_retained = if (n) 1 else NA_real_, tolerance = 0)
  list(W = rbind(Wi, Wb), Wi = Wi, Wb = Wb, interior = info(length(ii)),
       boundary = info(length(bb)), conditional_covariance = diag(length(bb)), ii = ii, bb = bb)
}

.jgp_H <- function(tt, obs) {
  n <- length(tt); H <- matrix(0, length(obs), n)
  left <- pmax(1L, pmin(n - 1L, findInterval(obs, tt)))
  alpha <- (obs - tt[left]) / (tt[left + 1L] - tt[left])
  H[cbind(seq_along(obs), left)] <- 1 - alpha
  H[cbind(seq_along(obs), left + 1L)] <- alpha
  H
}

.jgp_grid <- function(tt, control) {
  n <- length(tt); h <- mean(diff(tt))
  regular <- max(abs(diff(tt) - h)) <= 1e-8 * h
  # Use the smallest uniform observation-aligned grid meeting grid_min. A
  # power-of-two subdivision was only needed by the removed refinement loop.
  m <- if (regular) (n - 1L) * max(1, ceiling((control$grid_min - 1L) / (n - 1L))) + 1L
       else max(n, control$grid_min)
  if (m > control$grid_max) stop("grid_max cannot contain the fixed working grid.")
  seq(min(tt), max(tt), length.out = m)
}

.jgp_interpolation_error <- function(fit, tt, obs, H) {
  n <- length(obs); m <- length(tt)
  if (all(rowSums(H != 0) == 1L)) return(0)
  pred <- .jgp_gp_predict(fit, c(obs, tt))
  E <- cbind(diag(n), -H)
  bias <- as.vector(E %*% pred$mean)
  variance <- rowSums((E %*% pred$Sigma) * E)
  max(sqrt(pmax(0, bias^2 + variance))) / sqrt(fit$noise2)
}

# Scaled LM with explicit acceptance and a projected gradient for parameter
# bounds. Linear algebra uses QR on the augmented Jacobian, not its square.
.jgp_lm <- function(theta, evaluate, control, lower = rep(-Inf, length(theta)),
                    upper = rep(Inf, length(theta))) {
  implementation <- evaluate
  evaluations <- c(objective=0L,jacobian=0L)
  evaluate <- function(theta,jacobian=TRUE) {
    evaluations["objective"] <<- evaluations["objective"]+1L
    if (jacobian) evaluations["jacobian"] <<- evaluations["jacobian"]+1L
    implementation(theta,jacobian)
  }
  theta <- pmax(lower, pmin(upper, theta))
  current <- evaluate(theta, TRUE)
  if (any(!is.finite(current$r)) || any(!is.finite(current$J))) stop("Nonfinite initial LM residual/Jacobian.")
  value <- sum(current$r^2) / 2; damping <- control$damping
  history <- data.frame(iteration = 0L, objective = value, damping = damping,
                        accepted = TRUE, gradient = NA_real_)
  converged <- FALSE; reason <- "iteration limit"; iter <- 0L
  for (iter in seq_len(control$maxit)) {
    J <- current$J; gradient <- as.vector(crossprod(J, current$r))
    colscale <- pmax(sqrt(colSums(J^2)), 1e-8)
    projected <- gradient
    projected[theta <= lower & gradient > 0 | theta >= upper & gradient < 0] <- 0
    gn <- max(abs(projected) / colscale)
    if (gn <= control$gtol) { converged <- TRUE; reason <- "scaled gradient"; break }
    blocked <- (theta<=lower & gradient>0) | (theta>=upper & gradient<0)
    step <- if (any(blocked)) {
      free <- which(!blocked); A <- sweep(J[,free,drop=FALSE],2,colscale[free],"/")
      tryCatch({
        out <- numeric(length(theta))
        out[free] <- as.vector(qr.solve(rbind(A,diag(sqrt(damping),length(free))),
          c(-current$r,rep(0,length(free))),tol=1e-12))/colscale[free]
        out
      },error=function(e)NULL)
    } else if (!is.null(current$prior)) {
      # Exact Schur/Woodbury solve of the SAME damped normal equations. The
      # grid-state prior is diagonal in whitened coordinates, so dense QR over
      # all grid variables is unnecessary. No state modes are truncated here.
      tryCatch({
        meta <- current$prior; np <- meta$parameters
        ip <- seq_len(np); iz <- seq.int(np + 1L, length(theta))
        other <- J[-meta$rows, , drop = FALSE]
        P <- other[, ip, drop = FALSE]; Z <- other[, iz, drop = FALSE]
        diagonal <- meta$weight + damping * colscale[iz]^2
        if (ncol(Z) <= nrow(Z)) {
          # Factor the state system when it is smaller than the residual
          # system. It also avoids the subtractive Woodbury reconstruction.
          small <- chol(crossprod(Z)+diag(diagonal,length(diagonal)))
          state_solve <- function(B) .jgp_solve(small,B)
        } else {
          ZD <- sweep(Z, 2, diagonal, "/")
          small <- chol(diag(nrow(Z)) + tcrossprod(ZD, Z))
          state_solve <- function(B) {
            if (is.null(dim(B))) B <- matrix(B, ncol = 1L)
            DB <- B / diagonal
            DB - t(ZD) %*% .jgp_solve(small, Z %*% DB)
          }
        }
        cross <- crossprod(Z, P)
        hz <- state_solve(cbind(-gradient[iz], cross))
        reduced <- crossprod(P) + diag(damping * colscale[ip]^2, np) - crossprod(cross, hz[, -1L, drop = FALSE])
        rp <- -gradient[ip] - as.vector(crossprod(cross, hz[, 1L]))
        dp <- as.vector(.jgp_solve(chol((reduced + t(reduced)) / 2), rp))
        dz <- hz[, 1L] - hz[, -1L, drop = FALSE] %*% dp
        c(dp, as.vector(dz))
      }, error = function(e) NULL)
    } else {
      A <- sweep(J, 2, colscale, "/")
      tryCatch(as.vector(qr.solve(rbind(A, diag(sqrt(damping), length(theta))),
                     c(-current$r, rep(0, length(theta))), tol = 1e-12)) / colscale,
                     error = function(e) NULL)
    }
    if (is.null(step) && !is.null(current$prior) && !any(blocked)) {
      # Near rank deficiency, subtracting the Schur complement can lose
      # positive definiteness. Fall back to the unsquared augmented system.
      A <- sweep(J,2,colscale,"/")
      step <- tryCatch(as.vector(qr.solve(rbind(A,diag(sqrt(damping),length(theta))),
        c(-current$r,rep(0,length(theta))),tol=1e-12))/colscale,error=function(e)NULL)
    }
    if (is.null(step)) { damping <- damping * 10; next }
    candidate <- pmax(lower, pmin(upper, theta + step)); step <- candidate - theta
    pred <- -sum(gradient * step) - sum((J %*% step)^2) / 2
    trial <- tryCatch(evaluate(candidate, FALSE), error = function(e) NULL)
    newvalue <- if (is.null(trial) || any(!is.finite(trial$r))) Inf else sum(trial$r^2) / 2
    gain <- if (pred > 0) (value - newvalue) / pred else -Inf
    accepted <- is.finite(gain) && gain > 1e-4 && newvalue < value
    history <- rbind(history, data.frame(iteration = iter, objective = if (accepted) newvalue else value,
                                        damping = damping, accepted = accepted, gradient = gn))
    if (accepted) {
      small_f <- value - newvalue <= control$ftol * max(1, value)
      small_x <- sqrt(sum(step^2)) <= control$xtol * (control$xtol + sqrt(sum(theta^2)))
      theta <- candidate; value <- newvalue; current <- evaluate(theta, TRUE)
      if (any(!is.finite(current$J))) stop("Nonfinite LM Jacobian after an accepted step.")
      damping <- max(1e-12, damping * max(1 / 3, 1 - (2 * gain - 1)^3))
      if (small_f || small_x) {
        # Small accepted improvements are reported distinctly from stationarity.
        converged <- TRUE; reason <- if (small_x) "step tolerance" else "objective tolerance"; break
      }
    } else damping <- min(1e16, damping * 5)
    if (damping >= 1e16) { reason <- "damping limit"; break }
  }
  gradient <- as.vector(crossprod(current$J, current$r))
  projected <- gradient
  projected[theta <= lower & gradient > 0 | theta >= upper & gradient < 0] <- 0
  list(par = theta, objective = value, converged = converged, reason = reason,
       iterations = iter, evaluations = evaluations, history = history, gradient = gradient,
       scaled_gradient = max(abs(projected) / pmax(sqrt(colSums(current$J^2)), 1e-8)),
       residual = current$r, jacobian = current$J)
}

.jgp_initialize <- function(U, p0, weak, model, Sigma, control, lower, upper) {
  D <- model$D; J <- model$J; m <- nrow(U)
  ii <- unlist(lapply(seq_len(D), function(d) (d - 1L) * weak$K + seq_len(weak$ni)))
  p <- p0
  first <- .jgp_weak_eval(U, p, weak, model)
  if (model$affine) {
    G <- first$Jp[ii, , drop = FALSE]
    target <- as.vector(G %*% p) - first$r[ii]
    scales <- pmax(sqrt(colSums(G^2)), 1e-12)
    scaled <- sweep(G, 2, scales, "/")
    if (qr(scaled)$rank < J) stop("Interior parameter initialization is rank deficient.")
    p <- pmax(lower, pmin(upper, as.vector(qr.solve(scaled, target)) / scales))
  } else {
    ev <- function(p, jacobian) {
      v <- .jgp_weak_eval(U, p, weak, model, jacobian)
      list(r = v$r[ii], J = if (jacobian) v$Jp[ii, , drop = FALSE] else NULL)
    }
    p <- .jgp_lm(p, ev, control, lower, upper)$par
  }
  fixed_metric <- if (.jgp_ode_mode(control) == "gp_delta") NULL else
    .jgp_ode_metric(weak, D, control, rows = ii)
  for (iteration in seq_len(control$init_gls)) {
    base <- .jgp_weak_eval(U, p, weak, model)
    S <- if (is.null(fixed_metric)) .jgp_ode_metric(weak, D, control, base$Ju, Sigma, ii) else fixed_metric
    W <- if (.jgp_ode_mode(control) == "identity") diag(length(ii)) else
      .jgp_whitener(S, control$covariance_tol)$W
    ev <- function(p, jacobian) {
      v <- .jgp_weak_eval(U, p, weak, model, jacobian)
      list(r = as.vector(W %*% v$r[ii]), J = if (jacobian) W %*% v$Jp[ii, , drop = FALSE] else NULL)
    }
    next_p <- .jgp_lm(p, ev, control, lower, upper)$par
    change <- max(abs(next_p - p) / pmax(1, abs(p))); p <- next_p
    if (change < control$xtol) break
  }
  p
}

.jgp_state <- function(fits, tt, control) {
  pred <- lapply(fits, .jgp_gp_predict, tt = tt)
  chol <- lapply(pred, function(x) .jgp_chol(x$K, control$gp_jitter, "GP state covariance"))
  state <- list(mean = do.call(cbind, lapply(pred, `[[`, "mean")),
       mu = matrix(rep(vapply(fits, `[[`, numeric(1), "mean"), each = length(tt)), length(tt)),
       L = lapply(chol, function(c) t(c$R)),
       Sigma = Map(function(p, c) p$Sigma + diag(c$added, length(tt)), pred, chol),
       radius = do.call(cbind, lapply(pred, `[[`, "radius")),
       jitter = vapply(chol, `[[`, numeric(1), "added"),
       include_gp_prior = !isFALSE(control$include_gp_prior))
  use_prior <- state$include_gp_prior
  state$bounds <- attr(control,"state_bounds")
  state$mean <- .jgp_clip_state(state$mean,state$bounds)
  state$coordinates <- .jgp_bounded_coordinates(list(
    mu=if (use_prior) state$mu else state$mu*0,
    L=if (use_prior) state$L else rep(list(1),ncol(state$mu)),
    order=seq_len(ncol(state$mu))),state$bounds)
  state$prior_whitened <- use_prior && (is.null(state$bounds) || !any(state$bounds$bounded))
  state
}

.jgp_state_start <- function(state) {
  coordinates <- state$coordinates
  coordinates$pscale <- numeric()
  .jgp_pack_coordinates(numeric(),state$mean,coordinates)
}

.jgp_unpack <- function(theta,state,J) {
  coordinates <- state$coordinates
  m <- nrow(state$mu); D <- ncol(state$mu)
  z <- matrix(theta[-seq_len(J)],m,D); U <- coordinates$mu
  for (d in seq_len(D)) {
    L <- coordinates$L[[d]]
    U[,d] <- U[,d]+if (length(L)==1L) L*z[,d] else as.vector(L%*%z[,d])
  }
  list(p=theta[seq_len(J)],U=U,z=if (state$prior_whitened) z else NULL)
}

# Observed-coordinate preparation for the shared objective.
.jgp_objective <- function(Y,H,state,weak,model,weights,noise,lambda) {
  use_prior <- !isFALSE(state$include_gp_prior)
  coordinates <- state$coordinates; coordinates$pscale <- rep(1,model$J)
  .jgp_joint_objective(Y,H,weak,model,weights,noise,lambda,coordinates,state,
    if (use_prior) seq_len(model$D) else integer(),prior_whitened=state$prior_whitened)
}

# Grid quadrature only. Both call sites are on the grid path; gp_gauss returns
# before reaching them and gates resolution through weak$resolved instead.
.jgp_refinement <- function(state, p, weak, model, fits, weights, control) {
  if (identical(weak$integration, "gp_gauss"))
    stop("Grid refinement diagnostics do not apply to gp_gauss integration.")
  fine_tt <- seq(min(weak$tt), max(weak$tt), length.out = 2 * length(weak$tt) - 1L)
  fine_mean <- do.call(cbind, lapply(fits, function(f) .jgp_gp_predict(f, fine_tt)$mean))
  fine_weak <- .jgp_weak(fine_tt, weak$design, control, weak$maps)
  coarse <- .jgp_weak_eval(state$mean, p, weak, model, FALSE)$r
  fine <- .jgp_weak_eval(fine_mean, p, fine_weak, model, FALSE)$r
  change <- coarse - fine
  wrms <- function(W) if (!nrow(W)) 0 else sqrt(mean((W %*% change)^2))
  c(interior = wrms(weights$Wi), boundary = wrms(weights$Wb))
}

.jgp_solution_accuracy <- function(U, p, state, weak, model, fits, weights, H, control) {
  if (identical(weak$integration, "gp_gauss"))
    return(.jgp_gauss_accuracy(U, p, state, weak, model, fits, weights, H, control))
  fine_tt <- seq(min(weak$tt), max(weak$tt), length.out = 2 * length(weak$tt) - 1L)
  fine_U <- matrix(0, length(fine_tt), ncol(U))
  observation_error <- numeric(ncol(U))
  exact <- all(rowSums(H != 0) == 1L)
  for (d in seq_len(ncol(U))) {
    if (isFALSE(state$include_gp_prior)) {
      # A free grid state has the piecewise-linear representation used by H,
      # not a GP-conditioned extension that would reintroduce GP assumptions.
      fine_U[,d] <- .jgp_H(weak$tt,fine_tt) %*% U[,d]
      next
    }
    # GP conditional extension of the optimized grid state. Numerical jitter is
    # only on existing grid variables, not a new white-noise signal between them.
    alpha <- .jgp_solve(t(state$L[[d]]), U[, d] - state$mu[, d])
    fine_U[, d] <- fits[[d]]$mean + .jgp_gp_cross(fits[[d]], fine_tt, weak$tt) %*% alpha
    if (!exact) {
      observed <- fits[[d]]$mean + .jgp_gp_cross(fits[[d]], fits[[d]]$tt, weak$tt) %*% alpha
      observation_error[d] <- max(abs(observed - H %*% U[, d])) / sqrt(fits[[d]]$noise2)
    }
  }
  # A diagonal numerical stabilization lives only at the original grid nodes.
  # Its optimized contribution may be tiny in state units yet large after weak
  # whitening. Separate this from quadrature of the continuous GP extension;
  # otherwise increasing grid_min is misleading advice for a stabilization floor.
  fw <- .jgp_weak(fine_tt, weak$design, control, weak$maps)
  coarse_raw <- .jgp_weak_eval(U, p, weak, model, FALSE)$r
  smooth_U <- fine_U[seq(1L, length(fine_tt), by = 2L), , drop = FALSE]
  smooth_raw <- .jgp_weak_eval(smooth_U, p, weak, model, FALSE)$r
  smooth_fine_raw <- .jgp_weak_eval(fine_U, p, fw, model, FALSE)$r
  # Preserve exact values on the nested coarse grid, including its numerical
  # stabilization. This original check also sees any on-grid/extension mismatch;
  # the continuous-extension diagnostics above identify that contribution.
  fine_U[seq(1L, length(fine_tt), by = 2L), ] <- U
  delta <- coarse_raw - .jgp_weak_eval(fine_U, p, fw, model, FALSE)$r
  block_rms <- function(delta) {
    rms <- function(W) if (!nrow(W)) 0 else sqrt(mean((W %*% delta)^2))
    if (control$lambda > 0) c(interior = rms(weights$Wi), boundary = rms(weights$Wb))
    else c(interior = 0, boundary = 0)
  }
  errors <- block_rms(delta)
  list(interpolation = observation_error, weak = errors,
       smooth_quadrature = block_rms(smooth_raw - smooth_fine_raw),
       stabilization_weak = block_rms(coarse_raw - smooth_raw),
       stabilization_state = apply(abs(U - smooth_U), 2, max),
       weak_passed = max(errors) <= control$weak_grid_tol,
       passed = max(observation_error) <= control$grid_tol)
}

#' Joint weak-form GP estimation with observed or latent components
#'
#' @param f Function f(u,p,t), compatible with WENDy's symbolic machinery.
#' @param U Numeric observation matrix, one column per system component. The
#'   observed formulation requires all finite values. The latent formulation
#'   requires exactly one entirely NA column and finite observed columns.
#' @param tt Strictly increasing finite observation times.
#' @param p0 Optional initial parameter vector. Required if parameter indices
#'   cannot be inferred from f. Observed fits refine it by interior-only GLS.
#' @param noise_sd Known positive scalar or per-observed-component measurement SD, or
#'   NULL to estimate it by GP marginal likelihood.
#' @param control Named overrides from [wendygp_control()].
#' @param parameter_lower,parameter_upper Optional physical parameter bounds
#'   for both formulations. Supply vectors in parameter-index order (p[1],
#'   p[2], ...); names are labels, not used for matching. A scalar applies to
#'   every parameter. NULL means unbounded on that side. Each lower bound must
#'   be strictly smaller than its upper bound.
#' @param state_lower,state_upper Optional physical state bounds, scalar or
#'   one value per column of U in column order, including the latent column.
#'   These constrain fitted grid states, not the entire continuous GP curve.
#'   NULL means unbounded on that side.
#' @param lower,upper Compatibility aliases for parameter_lower and
#'   parameter_upper. Do not supply an alias and its corresponding explicit
#'   argument together.
#' @param formulation "auto" (default) selects "observed" for finite data or
#'   "latent" for exactly one entirely NA column with finite observed columns.
#'   Partial missingness and multiple latent columns are unsupported. Set
#'   "observed" or "latent" to require a particular formulation.
#' @param solver Optimizer: NULL selects LM for observed and nlminb for latent.
#'   Observed fits also accept "nlminb" and "lbfgsb"; latent fits accept those
#'   two scalar optimizers. All use analytic derivatives of the same objective.
#' @param starts Finite parameter multipliers for latent starts (default 1).
#'   Stationary completed runs are preferred, followed by finite feasible
#'   stage-3 estimates, then stage-1 fallbacks. Runs in the preferred group are
#'   ranked by objective value. latent_penalty=TRUE requires
#'   one start because different fitted priors are different objectives.
#' @param latent_initial Optional finite latent starting curve at observation
#'   times. Default: row mean of observations, falling back to their largest-RMS
#'   column if it cancels. Changing this curve does not change frozen ODE scales.
#' @param latent_penalty Add a frozen latent GP quadratic in stage 3 (default
#'   FALSE). Otherwise the fitted covariance is only a preconditioner.
#' @param maxit Optional iteration cap, overriding control$maxit. Defaults to
#'   500 for observed and 5000 per latent optimization phase.
#' @param gradient_tol Optional physical stationarity tolerance, overriding
#'   control$gtol. Defaults to 1e-6 for observed and 1e-4 for latent.
#' @param modes Latent compatibility argument: NULL or the full numerical span's
#'   size. Truncation is unsupported.
#' @param screen Latent compatibility argument; must be FALSE.
#' @details Fits nonstationary Matern 5/2 GPs by default, with profiled constant
#' means and signal variances, then freezes the hyperparameters. A bounded smooth
#' radius function supplies nonstationarity. Set control=list(kernel="bump")
#' for normalized bump-convolution covariance. Both kernels use bump weak tests.
#' By default, interior tests come
#' from a geometric multiscale pool independent of GP radii, with a test at every
#' admissible working-grid center per radius. Integration-accuracy screening
#' precedes adaptive SVD truncation at weak_design_info=0.95 (cumulative singular
#' values, not squares). The working grid and pool are constructed once; nested
#' quadrature checks never enlarge the optimization grid or restart the solve. Set
#' control=list(weak_radius_method="gp") for the legacy radius-transfer method.
#' Grid acceptance checks only observation-operator error. Weak quadrature
#' diagnostics are returned separately and never reject the fitted grid.
#' The default control=list(weak_integration="gp_gauss") evaluates continuous
#' GP weak integrals, their Jacobians and the full Gram consistently, with exact
#' endpoints and no Euler-Maclaurin substitution. The optimized grid is unchanged.
#' This path requires the GP prior and SVD or explicit weak tests;
#' its selected-space quadrature checks do not automatically change that space.
#' With control=list(include_gp_prior=FALSE), the objective contains only data
#' and lambda-weighted weak penalties. Optimization is directly over U;
#' gp_delta covariance stays frozen and sensitivity scores use no GP prior.
#' Free-state grid diagnostics use H's piecewise-linear representation, not a
#' GP-conditioned extension. A finite weak space need not identify every state
#' direction, so this option is not necessarily a proper posterior.
#' An interior-only frozen-state fit
#' initializes p. The default test-Gram ODE penalty uses RMS/span normalization
#' (ode_units="rms_span"); use ode_units="raw" for physical-unit Gram weighting.
#' control=list(ode_weighting="gp_delta") instead propagates the frozen posterior
#' grid-state covariance, including endpoint/EM channels. Interior and boundary
#' cross terms are retained by conditional whitening. ode_weighting="identity"
#' uses unweighted retained weak rows. The default does not propagate GP
#' covariance into the ODE penalty; the data likelihood and GP prior are unchanged.
#' problem$Omega holds the selected metric (a covariance only in gp_delta mode).
#' rho_initial and rho_final are mean squared metric-weighted interior residuals;
#' only gp_delta expresses them in GP-standardized units. The selected metric
#' and basis remain fixed during joint optimization. The likelihood and GP prior are
#' untempered; only the ODE penalty receives lambda, which defaults to 100.
#' Set control=list(weak_radius_method="posterior") to use the experimental
#' posterior/noise-floor weak-radius selector instead of the multiscale SVD pool.
#' Its profile and assumptions are returned in diagnostics$radius_selection.
#' An inadmissible automatic pool raises a wendygp_radius_error condition with
#' its radius_selection profile attached; no fallback radius is substituted.
#' The latent formulation uses the same fixed-covariance objective engine,
#' with likelihood rows only for observed components. Stage 1 optimizes physical
#' latent values with no latent prior. Stage 2 fits a GP to that curve. Stage 3
#' freezes this covariance and changes optimization coordinates while keeping
#' the stage-1 curve as its start. When latent_penalty=TRUE, it instead starts
#' from the smoothed curve and adds the frozen latent quadratic.
#' No stage jointly optimizes latent values and GP amplitude. With the quadratic
#' off, stages 1 and 3 have the same physical objective.
#' Latent fits require matern52, gp_gauss, observed GP priors, test_gram or
#' identity weighting, and full numerical test span. Their observed
#' extension is Lagrange and latent extension is GP. Default latent ODE scales
#' are initializer-based heuristics; use ode_component_scale for known units.
#' This is a plug-in generalized MAP estimate, not full Bayesian inference.
#' @return A jointgp object with phat, U_hat, tt, gp fits, initial estimates,
#' objective contributions, convergence diagnostics, and a reusable problem
#' object for derivative checks and controlled ablations. GP covariance
#' stabilization and projected weak-covariance ranks are reported explicitly.
#' formulation records the resolved formulation. converged reports
#' physical-coordinate stationarity. optimizer separately records the native
#' exit reason, convergence flag, iterations, and evaluation counts.
#' problem$evaluate(theta) returns residual r, Jacobian J, scalar value and
#' analytic gradient. With scalar=TRUE it avoids assembling the stacked
#' Jacobian; jacobian=FALSE requests values only. Latent fits additionally include
#' stage records in runs, observed/latent column indices, and class wendygp_latent.
#' A finite feasible latent stage-3 estimate is returned even without
#' stationarity, with a warning and converged=FALSE. Stage 1 is returned only
#' if no usable stage-3 estimate exists. If neither stage has a finite feasible
#' estimate, a wendygp_latent_numerical_error carries all runs.
#' solveWendyGPLatent is a compatibility wrapper.
#' diagnostics$stationarity includes a local joint Jacobian rank; neither that
#' rank nor stationarity guarantees global identifiability. On the Gauss path,
#' quadrature_passed and weak_grid_passed describe the returned estimate's
#' residual, Jacobian and Gram integration checks. Initial checks are retained
#' for inspection. extension and extension_passed report advisory GP conditional
#' uncertainty relative to fitted-state RMS; they do not gate numerical accuracy.
#' @examples
#' \dontrun{
#' tt <- seq(0, 8, length.out = 41)
#' set.seed(12)
#' U <- matrix(1 / (1 + 9 * exp(-tt)) + rnorm(length(tt), sd = 0.03))
#' fit <- solveWendyGP(function(u,p,t) c(p[1]*u[1]*(1-u[1])),
#'                     U, tt, p0 = 0.8)
#' fit$phat
#' }
#' @export
solveWendyGP <- function(f,U,tt,p0=NULL,noise_sd=NULL,control=NULL,
                         parameter_lower=NULL,parameter_upper=NULL,formulation=c("auto","observed","latent"),
                         solver=NULL,starts=1,latent_initial=NULL,latent_penalty=FALSE,
                         maxit=NULL,gradient_tol=NULL,modes=NULL,screen=FALSE,
                         state_lower=NULL,state_upper=NULL,lower=NULL,upper=NULL) {
  lower <- .jgp_parameter_bound_alias(parameter_lower,lower,"parameter_lower","lower")
  upper <- .jgp_parameter_bound_alias(parameter_upper,upper,"parameter_upper","upper")
  formulation <- match.arg(formulation)
  if (formulation=="auto") {
    if (is.vector(U) && is.numeric(U)) U <- matrix(U,ncol=1L)
    formulation <- "observed"
    if (is.matrix(U) && is.numeric(U) && nrow(U)>0L && ncol(U)>0L && anyNA(U)) {
      .jgl_split(U) # Reject partial missingness, multiple latent columns or no observations.
      formulation <- "latent"
    }
  }
  if (is.null(solver)) solver <- if (formulation=="observed") "lm" else "nlminb"
  solver <- match.arg(solver,c("lm","nlminb","lbfgsb"))
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control must be a list of wendygp_control overrides.")
  if (!is.null(maxit)) control$maxit <- maxit
  if (!is.null(gradient_tol)) control$gtol <- gradient_tol
  tt <- as.vector(tt)
  if (formulation=="latent") {
    if (solver=="lm") stop("The latent formulation requires solver='nlminb' or 'lbfgsb'.")
    if (is.null(control$maxit)) control$maxit <- 5000L
    if (is.null(control$gtol)) control$gtol <- 1e-4
    answer <- .jgp_latent_solve(f,U,tt,p0,noise_sd,control,starts,modes,screen,
      latent_penalty,maxit=control$maxit,gradient_tol=control$gtol,
      solver=solver,latent_initial=latent_initial,lower=lower,upper=upper,
      state_lower=state_lower,state_upper=state_upper)
  } else {
    if (!is.null(latent_initial) || !identical(latent_penalty,FALSE) ||
        !identical(as.numeric(starts),1) || !is.null(modes) || !identical(screen,FALSE))
      stop("starts, latent_initial, latent_penalty, modes and screen apply only to formulation='latent'.")
    answer <- .jgp_observed_solve(f,U,tt,p0,noise_sd,control,lower,upper,solver,state_lower,state_upper)
    answer$observed <- seq_len(ncol(answer$U_hat)); answer$latent <- integer()
  }
  answer$bounds <- list(parameter_lower=answer$problem$lower,parameter_upper=answer$problem$upper,
    state_lower=answer$problem$state_bounds$lower,state_upper=answer$problem$state_bounds$upper)
  answer$formulation <- formulation; answer$call <- match.call()
  answer
}

.jgp_observed_solve <- function(f, U, tt, p0 = NULL, noise_sd = NULL, control = NULL,
                         lower = NULL, upper = NULL, solver = "lm",
                         state_lower = NULL, state_upper = NULL) {
  control <- do.call(wendygp_control, if (is.null(control)) list() else control)
  attr(control,"solver") <- solver
  if (is.vector(U) && is.numeric(U)) U <- matrix(U, ncol = 1L)
  if (!is.matrix(U) || !is.numeric(U) || nrow(U) < 6L || ncol(U) < 1L)
    stop("U must be a numeric matrix with at least six observations per component.")
  if (any(!is.finite(U))) stop("The observed formulation requires finite observations in every column. Use formulation='latent' for one entirely NA column.")
  tt <- as.vector(tt)
  if (!is.numeric(tt) || length(tt) != nrow(U) || any(!is.finite(tt)) || any(diff(tt) <= 0))
    stop("tt must be finite, strictly increasing, and match U's rows.")
  D <- ncol(U); n <- nrow(U)
  attr(control,"state_bounds") <- .jgp_bound_pair(state_lower,state_upper,D,"state")
  J <- if (is.null(p0)) detect_n_params(f) else length(p0)
  if (J < 1L) stop("Supply p0 so the number of parameters is known.")
  if (is.null(p0)) p0 <- rep(1, J)
  if (!is.numeric(p0) || any(!is.finite(p0))) stop("p0 must be a finite numeric vector.")
  bounds <- .jgp_bound_pair(lower,upper,J,"parameter")
  lower <- bounds$lower; upper <- bounds$upper
  p0 <- pmax(lower,pmin(upper,p0))
  if (!is.null(noise_sd)) {
    if (!is.numeric(noise_sd) || !length(noise_sd) %in% c(1L, D) || any(!is.finite(noise_sd) | noise_sd <= 0))
      stop("noise_sd must be NULL or a positive scalar/per-component vector.")
    noise_sd <- rep_len(noise_sd, D)
  }
  if (control$verbose) message("Fitting ", D, " component GP(s).")
  fits <- lapply(seq_len(D), function(d) .jgp_gp_fit(tt, U[, d],
    if (is.null(noise_sd)) NULL else noise_sd[d], control))
  noise <- sqrt(vapply(fits, `[[`, numeric(1), "noise2"))
  if (control$weak_integration == "gp_gauss") {
    answer <- .jgp_gauss_solve(f, U, tt, p0, fits, control, lower, upper)
    answer$call <- match.call()
    return(answer)
  }
  model <- .jgp_model(f, D, J, control$em_order)
  if (is.null(control$weak_radii) && control$weak_radius_method %in% c("sensitivity", "svd")) {
    group <- .jgp_design_group(U, tt, p0, fits, model, control, lower, upper,
                              arms = control$weak_radius_method)
    answer <- group$fits[[control$weak_radius_method]]
    if (inherits(answer, "error")) stop(answer)
    if (!answer$diagnostics$design_passed)
      warning("No fully screened candidate space was available; using the reported diagnostic fallback.", call. = FALSE)
    if (!answer$diagnostics$grid_passed) {
      msg <- "Observation-operator error exceeds grid_tol on the fixed working grid; inspect interpolation diagnostics or explicitly choose a larger grid_min."
      if (control$grid_action == "error") stop(msg) else warning(msg, call. = FALSE)
    }
    answer$call <- match.call()
    return(answer)
  }
  grid <- .jgp_grid(tt, control)
  design <- .jgp_test_design(tt, fits, control)
  if (control$verbose) message("Preparing fixed grid with ", length(grid), " points.")
  state <- .jgp_state(fits, grid, control)
  H <- .jgp_H(grid, tt)
  interp_error <- vapply(fits, .jgp_interpolation_error, numeric(1), tt = grid, obs = tt, H = H)
  weak <- .jgp_weak(grid, design, control)
  pinit <- .jgp_initialize(state$mean, p0, weak, model, state$Sigma, control, lower, upper)
  base <- .jgp_weak_eval(state$mean, pinit, weak, model)
  Omega <- .jgp_ode_metric(weak, D, control, base$Ju, state$Sigma)
  metric <- .jgp_fit_metric(Omega, weak, U, tt, control)
  weights <- metric$weights
  residual_error <- if (control$lambda > 0) .jgp_refinement(state, pinit, weak, model, fits, weights, control)
                    else c(interior = 0, boundary = 0)
  ok <- max(interp_error) <= control$grid_tol
  grid_history <- list(list(points = length(grid), interpolation = interp_error,
                           weak = residual_error, weak_passed = max(residual_error) <= control$weak_grid_tol,
                           passed = ok))
  if (!ok) {
    msg <- paste0("Observation-operator error exceeds grid_tol at ", length(grid), " fixed grid points: interpolation/noise = ",
      signif(max(interp_error), 3),
      ". Choose grid_min explicitly or use grid_action='warn' for a diagnostic fit; no automatic refinement is performed.")
    if (control$grid_action == "error") stop(msg) else warning(msg, call. = FALSE)
  }
  m <- length(grid)
  theta0 <- c(pinit, .jgp_state_start(state))
  objective <- .jgp_objective(U, H, state, weak, model, weights, noise, control$lambda)
  if (control$verbose) message("Joint ", solver, ": ", length(theta0), " variables.")
  coordinates <- state$coordinates; coordinates$pscale <- rep(1,J)
  bounds <- .jgp_coordinate_bounds(coordinates,lower,upper,state$bounds)
  result <- .jgp_optimize(theta0,objective,control,bounds$lower,bounds$upper)
  point <- .jgp_unpack(result$par, state, J)
  final <- objective(result$par, TRUE)
  stationarity <- .jgp_fit_diagnostics(point$U,point$p,U,H,noise,weak,model,weights,control,
    state,if (state$include_gp_prior) seq_len(D) else integer(),lower=lower,upper=upper)
  final_accuracy <- .jgp_solution_accuracy(point$U, point$p, state, weak, model, fits, weights, H, control)
  if (!final_accuracy$passed) {
    msg <- paste0("Final-solution observation-operator error exceeds grid_tol at ", m,
      " points: interpolation/noise = ", signif(max(final_accuracy$interpolation), 3),
      ". Check the chosen grid_min and numerical GP stabilization; no automatic refinement is performed.")
    if (control$grid_action == "error") stop(msg) else warning(msg, call. = FALSE)
  }
  rho <- function(r) sum((weights$Wi %*% r)^2) / weights$interior$rank
  structure(list(phat = point$p, U_hat = point$U, U_obs_hat = H %*% point$U,
    tt = grid, tt_obs = tt, Y = U, gp = fits, noise_sd = noise,
    initial = list(p = pinit, U = state$mean, theta = theta0),
    lambda = control$lambda, ode_weighting = control$ode_weighting,
    ode_units = metric$scaling$units,
    include_gp_prior = state$include_gp_prior,
    objective = result$objective,
    contributions = final$contributions, converged = stationarity$stationary,
    optimizer = list(method=result$method,converged=result$converged,reason=result$reason,
      iterations=result$iterations,evaluations=result$evaluations),
    convergence_reason = result$reason, iterations = result$iterations,
    diagnostics = list(rho_initial = rho(base$r), rho_final = rho(final$raw),
      ode_weighting = control$ode_weighting,
      ode_scaling = metric$scaling,
      include_gp_prior = state$include_gp_prior,
      state_coordinates = if (state$prior_whitened) "gp_whitened" else "physical_or_mixed",
      gp_uncertainty_propagated = control$ode_weighting == "gp_delta",
      rho_is_gp_standardized = control$ode_weighting == "gp_delta", ode_metric_frozen = TRUE,
      radius_selection = design$radius_selection,
      stationarity = stationarity,
      scaled_gradient = stationarity$scaled_gradient, gradient = stationarity$gradient,
      interior_rank = weights$interior$rank, boundary_rank = weights$boundary$rank,
      interior_discarded = weights$interior$discarded, boundary_discarded = weights$boundary$discarded,
      gp_jitter = state$jitter, gp_converged = vapply(fits, function(g) g$convergence == 0, logical(1)),
      grid_passed = ok && final_accuracy$passed, grid_history = grid_history,
      grid_fixed = TRUE, grid_refinements = 0L,
      grid_acceptance = "observation_operator",
      weak_grid_passed = max(residual_error) <= control$weak_grid_tol && final_accuracy$weak_passed,
      final_grid = final_accuracy, history = result$history,
      covariance_frozen = TRUE, hyperparameters_frozen = TRUE),
    problem = list(evaluate = objective, theta = result$par, state = state, H = H,
      weak = weak, model = model, Omega = metric$Omega, Omega_raw = Omega,
      weights = weights, control = control,
      lower = lower, upper = upper,state_bounds=state$bounds), call = match.call()), class = "jointgp")
}

#' @export
print.jointgp <- function(x, ...) {
  cat("WENDyGP (", if (is.null(x$formulation)) "observed" else x$formulation,
      "; ", x$gp[[1]]$control$kernel, " covariance)\n", sep = "")
  cat("Parameters:", format(x$phat, digits = 6), "\n")
  cat("Objective:", format(x$objective, digits = 6), " lambda:", x$lambda, "\n")
  cat("ODE weighting:", .jgp_ode_mode(x$problem$control), "\n")
  cat("ODE units:", if (is.null(x$ode_units)) "raw" else x$ode_units, "\n")
  if (identical(x$weak_integration,"gp_gauss"))
    cat("Weak integration: continuous GP / Gauss (fixed optimized grid)\n")
  cat("GP prior:", if (isFALSE(x$include_gp_prior)) "disabled" else "enabled", "\n")
  cat("Converged:", x$converged, "\n")
  if (!is.null(x$diagnostics$stationarity))
    cat("Physical stationarity:", x$diagnostics$stationarity$stationary, "\n")
  cat("Optimizer exit:", if (is.null(x$optimizer$reason)) x$convergence_reason else x$optimizer$reason, "\n")
  if (!x$converged) cat("Fit status:", x$convergence_reason, "\n")
  cat("Observation-operator accuracy passed:", x$diagnostics$grid_passed, "\n")
  invisible(x)
}
