# WENDyGP interface: solveWendyGP (input checks and dispatch to the fits in
# joint_gp_fit.R), wendygp_control, and printing. Stacked vectors are
# component-major.

#' Joint weak-form GP estimation with observed or latent components
#'
#' @param f Function f(u,p,t), compatible with WENDy's symbolic machinery.
#' @param U Numeric observation matrix, one column per system component. The
#'   observed formulation requires all finite values. The latent formulation
#'   requires exactly one entirely NA column and finite observed columns.
#' @param tt Strictly increasing finite observation times.
#' @param p0 Optional initial parameter vector. Required if the parameter count
#'   cannot be inferred from f.
#' @param noise_sd Known positive scalar or per-observed-component measurement
#'   SD, or NULL to estimate it by GP marginal likelihood.
#' @param control Named overrides from [wendygp_control()].
#' @param parameter_lower,parameter_upper Optional parameter bounds in
#'   parameter-index order; a scalar applies to every parameter. NULL means
#'   unbounded on that side.
#' @param formulation "auto" (default) selects "latent" when U has NA values
#'   and "observed" otherwise.
#' @param latent_initial Optional finite latent starting curve at observation
#'   times. Default: row mean of the observed columns (largest-RMS column if it
#'   cancels). It does not change the frozen ODE scales or latent GP.
#' @param maxit Optional iteration cap, overriding control$maxit. Defaults to
#'   500 for observed and 5000 per latent optimization.
#' @param gradient_tol Optional optimizer gradient tolerance, overriding
#'   control$gtol. Defaults to 1e-6 for observed and 1e-4 for latent.
#' @param state_lower,state_upper Optional state bounds, scalar or one value
#'   per column of U (including the latent column). They constrain the fitted
#'   grid states.
#' @details Observed components get nonstationary Matern 5/2 GPs whose
#' hyperparameters are then frozen; they set the noise level, starting state
#' and optimizer coordinates. The objective is the observation likelihood plus
#' lambda times the test-Gram-weighted weak ODE residual, with the ODE weight
#' scaled by n times the mean squared signal-to-noise ratio so it keeps pace
#' with the data. There is no GP prior term. Tests are multiscale bumps
#' orthonormalized on the grid. The residual is integrated by Gauss-Legendre
#' quadrature on each grid interval, with observed states extended to the nodes
#' by Lagrange interpolation. Parameters start from a frozen-state fit; states
#' and parameters are then fitted jointly by projected Levenberg-Marquardt.
#'
#' In the latent formulation the observed components keep GP prior rows, lambda
#' is not scaled by the data, and the hidden component has no data or prior rows. It
#' is interpolated by a fixed stationary Matern 5/2 with length scale
#' 1/gamma_t, where gamma_t = max|phi'| / (largest test radius) is the WSINDy
#' scale (arXiv:2007.02848, eq. 4.7). Default latent units assume comparable
#' component units; set ode_component_scale for known units. The estimate is a
#' plug-in MAP, not full Bayesian inference.
#' @return A jointgp object: phat; U_hat on the grid tt; tt_obs and Y; the
#'   observed and latent column indices; the observed GP fits gp and noise_sd;
#'   objective and contributions; converged (the optimizer met a tolerance) and
#'   convergence_reason; optimizer; bounds; formulation; and problem, whose
#'   evaluate(theta, jacobian = TRUE) returns the residual, Jacobian, value and
#'   gradient. Latent fits add class wendygp_latent. A nonconverged fit warns.
#' @examples
#' \dontrun{
#' tt <- seq(0, 8, length.out = 41)
#' set.seed(12)
#' U <- matrix(1 / (1 + 9 * exp(-tt)) + rnorm(length(tt), sd = 0.03))
#' fit <- solveWendyGP(function(u, p, t) c(p[1] * u[1] * (1 - u[1])),
#'   U, tt,
#'   p0 = 0.8
#' )
#' fit$phat
#' }
#' @export
solveWendyGP <- function(f, U, tt, p0 = NULL, noise_sd = NULL, control = NULL,
                         parameter_lower = NULL, parameter_upper = NULL, formulation = c("auto", "observed", "latent"),
                         latent_initial = NULL,
                         maxit = NULL, gradient_tol = NULL,
                         state_lower = NULL, state_upper = NULL) {
  if (is.vector(U) && is.numeric(U)) U <- matrix(U, ncol = 1L)
  if (!is.matrix(U) || !is.numeric(U) || nrow(U) < 6L || ncol(U) < 1L) {
    stop("U must be a numeric matrix with at least six rows.")
  }
  formulation <- match.arg(formulation)
  if (formulation == "auto") formulation <- if (anyNA(U)) "latent" else "observed"
  latent_fit <- formulation == "latent"
  if (latent_fit) {
    components <- identify_latent_component(U)
  } else {
    if (any(!is.finite(U))) {
      stop(
        "The observed formulation requires finite observations in every column. ",
        "Use formulation='latent' for one entirely NA column."
      )
    }
    components <- list(observed = seq_len(ncol(U)), latent = integer())
  }
  tt <- as.vector(tt)
  if (!is.numeric(tt) || length(tt) != nrow(U) || any(!is.finite(tt)) || any(diff(tt) <= 0)) {
    stop("tt must be finite, strictly increasing, and match U's rows.")
  }

  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control must be a list of wendygp_control overrides.")
  if (!is.null(maxit)) control$maxit <- maxit
  if (!is.null(gradient_tol)) control$gtol <- gradient_tol
  if (latent_fit) {
    if (is.null(control$maxit)) control$maxit <- 5000L
    if (is.null(control$gtol)) control$gtol <- 1e-4
    if (is.null(control$lambda)) control$lambda <- 100
  }
  control <- do.call(wendygp_control, control)
  n_components <- ncol(U)
  state_bounds <- validate_joint_gp_bounds(state_lower, state_upper, n_components, "state")
  attr(control, "state_bounds") <- state_bounds

  n_parameters <- if (is.null(p0)) detect_n_params(f) else length(p0)
  if (n_parameters < 1L) stop("Supply p0 so the number of parameters is known.")
  if (is.null(p0)) p0 <- rep(1, n_parameters)
  if (!is.numeric(p0) || any(!is.finite(p0))) stop("p0 must be a finite numeric vector.")
  parameter_bounds <- validate_joint_gp_bounds(parameter_lower, parameter_upper, n_parameters, "parameter")
  lower <- parameter_bounds$lower
  upper <- parameter_bounds$upper
  p0 <- pmax(lower, pmin(upper, p0))
  n_observed <- length(components$observed)
  if (!is.null(noise_sd)) {
    if (!is.numeric(noise_sd) || !length(noise_sd) %in% c(1L, n_observed) || any(!is.finite(noise_sd) | noise_sd <= 0)) {
      stop("noise_sd must be NULL or positive, scalar or one value per observed component.")
    }
    noise_sd <- rep_len(noise_sd, n_observed)
  }

  if (latent_fit) {
    if (!is.null(latent_initial) && (!is.numeric(latent_initial) ||
      length(latent_initial) != nrow(U) || any(!is.finite(latent_initial)))) {
      stop("latent_initial must be finite with one value per observation time.")
    }
    answer <- fit_latent_joint_gp(f, U, tt, p0, noise_sd, control, latent_initial, lower, upper)
  } else {
    if (!is.null(latent_initial)) stop("latent_initial applies only to formulation='latent'.")
    answer <- fit_observed_joint_gp(f, U, tt, p0, noise_sd, control, lower, upper)
  }
  answer$bounds <- list(
    parameter_lower = lower, parameter_upper = upper,
    state_lower = state_bounds$lower, state_upper = state_bounds$upper
  )
  answer$formulation <- formulation
  answer$call <- match.call()
  answer
}

#' Controls for WENDyGP
#'
#' @param ... Named overrides. See Details for the available controls.
#' @details Defaults are returned by calling this function without arguments.
#' Radii are in physical time units; NULL selects data-based defaults.
#'
#' GP: nonstationary Matern 5/2 per observed component. gp_nonstationary_terms
#' cosine terms (0 for stationary) shape its log radius, penalized by
#' gp_radius_penalty and bounded by gp_radius_bounds. gp_restarts and gp_maxit
#' control the marginal-likelihood fit; gp_jitter is relative stabilization.
#'
#' Grid: uniform, with at least grid_min points (default 64) and two intervals
#' per observation interval; regular observation times stay on it. grid_max
#' (default 1000) caps its size.
#'
#' Tests: bumps (shape bump_eta) at every admissible grid point for each radius
#' in test_radii (NULL halves 0.4 * span down to one grid spacing),
#' orthonormalized on the grid keeping directions above basis_tol.
#' include_boundary adds boundary_count tests per end with radii boundary_radii
#' (default: the largest radius). covariance_tol is the whitening rank tolerance.
#'
#' Integration: quad_order Gauss-Legendre points per grid interval (default 8).
#' parameter_scale fixes the latent optimizer's parameter scales (default
#' abs(p0), 1 for zeros).
#'
#' ODE penalty: test Gram metric with component d scaled by T / a_d^2 (time
#' span over squared RMS, or ode_component_scale). In the observed formulation
#' it is further multiplied by n * mean(rms^2 / noise^2), so lambda (default
#' 0.3) is relative to the data likelihood. The latent formulation uses lambda
#' unscaled (default 100).
#'
#' Optimization: projected Levenberg-Marquardt, capped at maxit iterations. It
#' converges when the column-scaled projected gradient is below gtol or the
#' objective or scaled step changes by less than ftol or xtol. init_gls sets
#' the frozen-state passes of the parameter pilot. verbose prints progress.
#' @return A named list of controls.
#' @export
wendygp_control <- function(...) {
  defaults <- list(
    bump_eta = 9,
    gp_nonstationary_terms = 2L,
    gp_radius_bounds = NULL,
    gp_radius_penalty = 1,
    gp_restarts = 1L,
    gp_maxit = 200L,
    gp_jitter = 1e-10,
    grid_min = 64L,
    grid_max = 1000L,
    test_radii = NULL,
    include_boundary = TRUE,
    boundary_count = 1L,
    boundary_radii = NULL,
    basis_tol = 1e-7,
    covariance_tol = 1e-9,
    quad_order = 8L,
    parameter_scale = NULL,
    ode_component_scale = NULL,
    lambda = 0.3,
    init_gls = 5L,
    maxit = 500L,
    ftol = 1e-9,
    xtol = 1e-8,
    gtol = 1e-6,
    verbose = FALSE
  )
  args <- list(...)
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))) ||
    anyDuplicated(names(args)))) {
    stop("Controls must have unique names.")
  }
  unknown <- setdiff(names(args), names(defaults))
  if (length(unknown)) stop("Unknown WENDyGP control(s): ", paste(unknown, collapse = ", "))
  for (name in names(args)) {
    defaults[name] <- args[name]
  }
  control <- defaults
  positive <- c(
    "bump_eta", "gp_restarts", "gp_maxit", "gp_jitter", "grid_min", "grid_max",
    "boundary_count", "basis_tol", "covariance_tol", "quad_order", "maxit", "ftol", "xtol", "gtol"
  )
  for (name in positive) {
    if (!is.numeric(control[[name]]) || length(control[[name]]) != 1L ||
      !is.finite(control[[name]]) || control[[name]] <= 0) {
      stop(name, " must be one finite positive number.")
    }
  }
  for (name in c("gp_nonstationary_terms", "gp_radius_penalty", "init_gls", "lambda")) {
    if (!is.numeric(control[[name]]) || length(control[[name]]) != 1L ||
      !is.finite(control[[name]]) || control[[name]] < 0) {
      stop(name, " must be nonnegative.")
    }
  }
  integer_controls <- c(
    "gp_nonstationary_terms", "gp_restarts", "gp_maxit", "grid_min", "grid_max",
    "boundary_count", "init_gls", "maxit", "quad_order"
  )
  for (name in integer_controls) {
    if (control[[name]] != as.integer(control[[name]])) stop(name, " must be an integer.")
  }
  if (control$grid_min < 9 || control$grid_max < control$grid_min) stop("Require 9 <= grid_min <= grid_max.")
  if (control$quad_order < 2L || control$quad_order > 32L) stop("quad_order must be an integer from 2 to 32.")
  for (name in c("include_boundary", "verbose")) {
    if (!is.logical(control[[name]]) || length(control[[name]]) != 1L || is.na(control[[name]])) stop(name, " must be logical.")
  }
  for (name in c("gp_radius_bounds", "test_radii", "boundary_radii", "parameter_scale", "ode_component_scale")) {
    if (!is.null(control[[name]]) && (!is.numeric(control[[name]]) || !length(control[[name]]) ||
      any(!is.finite(control[[name]]) | control[[name]] <= 0))) {
      stop(name, " must be positive.")
    }
  }
  if (!is.null(control$gp_radius_bounds) && (length(control$gp_radius_bounds) != 2L ||
    diff(control$gp_radius_bounds) <= 0)) {
    stop("Invalid gp_radius_bounds.")
  }
  control
}

#' @export
print.jointgp <- function(x, ...) {
  control <- x$problem$control
  cat("WENDyGP (", x$formulation, ")\n", sep = "")
  cat("Parameters:", format(x$phat, digits = 6), "\n")
  cat("Objective:", format(x$objective, digits = 6), " lambda:", control$lambda, "\n")
  cat("Converged:", x$converged, "\n")
  cat("Optimizer exit:", x$optimizer$reason, "\n")
  if (!x$converged) cat("Fit status:", x$convergence_reason, "\n")
  invisible(x)
}
