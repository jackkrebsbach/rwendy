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
#' Final quadrature checks determine quadrature_passed and weak_grid_passed.
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
#' weak_design_radii supplies physical radii or NULL halves 0.4 times the span
#' down to one working-grid spacing. weak_design_quad_tol bounds relative quadrature
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
#' Trajectory and ODE parameter fits use projected Levenberg-Marquardt with
#' analytic Jacobians. LM adapts its damping internally; there is no damping
#' setting. GP marginal likelihood uses stats::nlminb with analytic gradients.
#' maxit sets the joint iteration budget; gp_maxit sets the GP-fitting budget.
#' ftol and xtol set objective and step tolerances. gtol tests the
#' Jacobian-column-scaled projected gradient during LM and physical stationarity
#' separately at the final fit.
#' init_gls controls frozen-state GLS iterations; xtol also stops these when
#' parameters stabilize. bump_eta fixes the bump shape. verbose prints progress.
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
    "bump_eta", "gp_restarts", "gp_maxit", "gp_quad_per_radius",
    "gp_quad_tol", "gp_jitter", "grid_min", "grid_max", "grid_tol",
    "weak_grid_tol", "weak_count", "bl_count", "basis_tol",
    "weak_radius_grid", "weak_radius_ratio", "weak_radius_tolerance",
    "weak_design_quad_tol", "weak_quad_order", "weak_quad_gram_tol",
    "covariance_tol", "maxit", "ftol", "xtol", "gtol"
  )
  for (name in positive) {
    if (!is.numeric(control[[name]]) || length(control[[name]]) != 1L ||
      !is.finite(control[[name]]) || control[[name]] <= 0) {
      stop(name, " must be one finite positive number.")
    }
  }
  for (name in c(
    "gp_nonstationary_terms", "gp_quad_max_refine", "gp_radius_penalty",
    "init_gls", "lambda", "weak_design_coverage"
  )) {
    if (!is.numeric(control[[name]]) || length(control[[name]]) != 1L ||
      !is.finite(control[[name]]) || control[[name]] < 0) {
      stop(name, " must be nonnegative.")
    }
  }
  integer_controls <- c(
    "gp_nonstationary_terms", "gp_quad_max_refine", "gp_restarts", "gp_maxit",
    "grid_min", "grid_max", "weak_count", "bl_count", "init_gls", "maxit",
    "weak_radius_grid", "weak_design_coverage", "weak_quad_order"
  )
  for (name in integer_controls) {
    if (control[[name]] != as.integer(control[[name]])) stop(name, " must be an integer.")
  }
  for (name in c("weak_design_budget", "weak_design_centers")) {
    if (!is.null(control[[name]]) && (!is.numeric(control[[name]]) || length(control[[name]]) != 1L ||
      !is.finite(control[[name]]) || control[[name]] <= 0 || control[[name]] > .Machine$integer.max ||
      control[[name]] != as.integer(control[[name]]))) {
      stop(name, " must be NULL or a positive integer.")
    }
  }
  if (!is.numeric(control$weak_design_info) || length(control$weak_design_info) != 1L ||
    !is.finite(control$weak_design_info) || control$weak_design_info <= 0 || control$weak_design_info > 1) {
    stop("weak_design_info must be in (0, 1].")
  }
  if (control$grid_min < 9 || control$grid_max < control$grid_min) stop("Require 9 <= grid_min <= grid_max.")
  if (control$weak_radius_grid < 9L) stop("weak_radius_grid must be at least 9.")
  if (control$weak_radius_ratio <= 1) stop("weak_radius_ratio must exceed 1.")
  if (!is.character(control$weak_radius_method) || length(control$weak_radius_method) != 1L ||
    is.na(control$weak_radius_method) || !control$weak_radius_method %in% c("gp", "posterior", "sensitivity", "svd")) {
    stop("weak_radius_method must be gp, posterior, sensitivity or svd.")
  }
  if (!is.null(control$weak_design_budget) && control$weak_design_coverage > control$weak_design_budget) {
    stop("weak_design_coverage cannot exceed weak_design_budget.")
  }
  if (!is.character(control$ode_weighting) || length(control$ode_weighting) != 1L ||
    is.na(control$ode_weighting) || !control$ode_weighting %in% c("gp_delta", "test_gram", "identity")) {
    stop("ode_weighting must be gp_delta, test_gram or identity.")
  }
  if (!is.character(control$ode_units) || length(control$ode_units) != 1L || is.na(control$ode_units) ||
    !control$ode_units %in% c("raw", "rms_span")) {
    stop("ode_units must be raw or rms_span.")
  }
  if (control$ode_weighting != "test_gram") {
    if ("ode_units" %in% names(args) && control$ode_units != "raw") {
      stop("ode_units='rms_span' requires ode_weighting='test_gram'.")
    }
    control$ode_units <- "raw"
  }
  if (!control$kernel %in% c("bump", "matern52")) stop("kernel must be bump or matern52.")
  if (!control$grid_action %in% c("error", "warn")) stop("grid_action must be error or warn.")
  if (length(control$em_order) != 1L || !control$em_order %in% c(0L, 2L, 4L)) stop("em_order must be 0, 2 or 4.")
  for (name in c("include_bl", "include_gp_prior", "verbose")) {
    if (!is.logical(control[[name]]) || length(control[[name]]) != 1L || is.na(control[[name]])) stop(name, " must be logical.")
  }
  if (!is.character(control$weak_extension) || !length(control$weak_extension) ||
    anyNA(control$weak_extension) || !all(control$weak_extension %in% c("lagrange", "gp"))) {
    stop("weak_extension must be lagrange or gp, one entry or one per component.")
  }
  if (!is.character(control$weak_integration) || length(control$weak_integration) != 1L ||
    is.na(control$weak_integration) || !control$weak_integration %in% c("grid", "gp_gauss")) {
    stop("weak_integration must be grid or gp_gauss.")
  }
  if (control$weak_quad_order < 2L || control$weak_quad_order > 32L) {
    stop("weak_quad_order must be an integer from 2 to 32.")
  }
  # The free-state and non-SVD formulations use grid quadrature.
  if (control$weak_integration == "gp_gauss" && (!control$include_gp_prior ||
    (is.null(control$weak_radii) && control$weak_radius_method != "svd"))) {
    control$weak_integration <- "grid"
    warning("weak_integration coerced to grid: gp_gauss supports neither a free ",
      "state (include_gp_prior=FALSE) nor weak_radius_method other than ",
      "svd without explicit weak_radii.",
      call. = FALSE
    )
  }
  for (name in c(
    "weak_radii", "weak_factors", "gp_radius_bounds", "bl_radii",
    "weak_design_radii", "weak_design_scale", "ode_component_scale"
  )) {
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
