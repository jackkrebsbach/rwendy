# WENDyGP public entry point and observed preparation. Both formulations use
# joint_gp_objective.R; latent staging lives in joint_gp_latent.R.
# Physical vectorizations are component-major (R column order).

build_observation_matrix <- function(grid, observation_times) {
  n_grid <- length(grid)
  observation_matrix <- matrix(0, length(observation_times), n_grid)
  left <- pmax(1L, pmin(n_grid - 1L, findInterval(observation_times, grid)))
  fraction <- (observation_times - grid[left]) / (grid[left + 1L] - grid[left])
  observation_matrix[cbind(seq_along(observation_times), left)] <- 1 - fraction
  observation_matrix[cbind(seq_along(observation_times), left + 1L)] <- fraction
  observation_matrix
}

build_joint_gp_grid <- function(tt, control) {
  n_observations <- length(tt)
  spacing <- mean(diff(tt))
  regular <- max(abs(diff(tt) - spacing)) <= 1e-8 * spacing
  # Use the smallest uniform observation-aligned grid meeting grid_min.
  n_grid <- if (regular) {
    (n_observations - 1L) * max(1, ceiling((control$grid_min - 1L) / (n_observations - 1L))) + 1L
  } else {
    max(n_observations, control$grid_min)
  }
  if (n_grid > control$grid_max) stop("grid_max cannot contain the fixed working grid.")
  seq(min(tt), max(tt), length.out = n_grid)
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
#' @param starts Finite parameter multipliers for latent starts (default 1).
#'   Stationary completed runs are preferred, followed by finite feasible
#'   stage-3 estimates, then stage-1 fallbacks. Runs in the preferred group are
#'   ranked by objective value. latent_penalty=TRUE requires
#'   one start because different fitted priors are different objectives.
#' @param latent_initial Optional finite latent starting curve at observation
#'   times. Default: row mean of observations, falling back to their largest-RMS
#'   column if it cancels. Changing this curve does not change frozen ODE scales.
#' @param latent_penalty Add a frozen latent GP quadratic in stage 3 (default
#'   FALSE). Otherwise the fitted covariance only changes optimization coordinates
#'   for unbounded states; bounded states use component scaling.
#' @param maxit Optional iteration cap, overriding control$maxit. Defaults to
#'   500 for observed and 5000 per latent optimization phase.
#' @param gradient_tol Optional physical stationarity tolerance, overriding
#'   control$gtol. Defaults to 1e-6 for observed and 1e-4 for latent.
#' @details Trajectory and ODE parameter fits use projected Levenberg-Marquardt
#' with analytic Jacobians and physical bounds. GP marginal likelihood uses
#' stats::nlminb with analytic gradients. Fits nonstationary Matern 5/2 GPs by default, with profiled constant
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
#' latent values with no latent prior. Stage 2 fits the latent GP mean, variance,
#' radius coefficients, and nugget to that curve. The nugget absorbs roughness
#' in the fitted curve and is excluded from the latent process covariance.
#' Stage 3 freezes the GP and changes optimization coordinates for unbounded
#' states while keeping the stage-1 curve as its start. Bounded states continue
#' to use component scaling. The weak interpolation and ODE weights stay fixed
#' across stages. When latent_penalty=TRUE, stage 3 instead starts
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
#' physical-coordinate stationarity. optimizer separately records the solver's
#' exit reason, convergence flag, iterations, and evaluation counts.
#' problem$evaluate(theta) returns residual r, Jacobian J, scalar value and
#' analytic gradient. jacobian=FALSE requests values only. Latent fits additionally include
#' stage records in runs, the separately fitted GP in latent_prior,
#' observed/latent column indices, and class wendygp_latent.
#' A finite feasible latent stage-3 estimate is returned even without
#' stationarity, with a warning and converged=FALSE. Stage 1 is returned only
#' if no usable stage-3 estimate exists. If neither stage has a finite feasible
#' estimate, or another numerical error occurs, the solve stops.
#' solveWendyGPLatent is a compatibility wrapper.
#' diagnostics$stationarity includes a local joint Jacobian rank; neither that
#' rank nor stationarity guarantees global identifiability. On the Gauss path,
#' quadrature_passed and weak_grid_passed describe the returned estimate's
#' residual, Jacobian and Gram integration checks. Initial checks are retained
#' for inspection.
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
                         starts = 1, latent_initial = NULL, latent_penalty = FALSE,
                         maxit = NULL, gradient_tol = NULL,
                         state_lower = NULL, state_upper = NULL, lower = NULL, upper = NULL) {
  lower <- resolve_parameter_bound_alias(parameter_lower, lower, "parameter_lower", "lower")
  upper <- resolve_parameter_bound_alias(parameter_upper, upper, "parameter_upper", "upper")
  formulation <- match.arg(formulation)
  if (formulation == "auto") {
    if (is.vector(U) && is.numeric(U)) U <- matrix(U, ncol = 1L)
    formulation <- "observed"
    if (is.matrix(U) && is.numeric(U) && nrow(U) > 0L && ncol(U) > 0L && anyNA(U)) {
      identify_latent_component(U) # Reject partial missingness, multiple latent columns or no observations.
      formulation <- "latent"
    }
  }
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control must be a list of wendygp_control overrides.")
  if (!is.null(maxit)) control$maxit <- maxit
  if (!is.null(gradient_tol)) control$gtol <- gradient_tol
  tt <- as.vector(tt)
  if (formulation == "latent") {
    if (is.null(control$maxit)) control$maxit <- 5000L
    if (is.null(control$gtol)) control$gtol <- 1e-4
    answer <- solve_latent_joint_gp(f, U, tt, p0, noise_sd, control, starts,
      latent_penalty,
      maxit = control$maxit, gradient_tol = control$gtol,
      latent_initial = latent_initial, lower = lower, upper = upper,
      state_lower = state_lower, state_upper = state_upper
    )
  } else {
    if (!is.null(latent_initial) || !identical(latent_penalty, FALSE) ||
      !identical(as.numeric(starts), 1)) {
      stop("starts, latent_initial and latent_penalty apply only to formulation='latent'.")
    }
    answer <- solve_observed_joint_gp(f, U, tt, p0, noise_sd, control, lower, upper, state_lower, state_upper)
    answer$observed <- seq_len(ncol(answer$U_hat))
    answer$latent <- integer()
  }
  answer$bounds <- list(
    parameter_lower = answer$problem$lower, parameter_upper = answer$problem$upper,
    state_lower = answer$problem$state_bounds$lower, state_upper = answer$problem$state_bounds$upper
  )
  answer$formulation <- formulation
  answer$call <- match.call()
  answer
}

solve_observed_joint_gp <- function(f, U, tt, p0 = NULL, noise_sd = NULL, control = NULL,
                                    lower = NULL, upper = NULL,
                                    state_lower = NULL, state_upper = NULL) {
  control <- do.call(wendygp_control, if (is.null(control)) list() else control)
  if (is.vector(U) && is.numeric(U)) U <- matrix(U, ncol = 1L)
  if (!is.matrix(U) || !is.numeric(U) || nrow(U) < 6L || ncol(U) < 1L) {
    stop("U must be a numeric matrix with at least six observations per component.")
  }
  if (any(!is.finite(U))) {
    stop(
      "The observed formulation requires finite observations in every column. ",
      "Use formulation='latent' for one entirely NA column."
    )
  }
  tt <- as.vector(tt)
  if (!is.numeric(tt) || length(tt) != nrow(U) || any(!is.finite(tt)) || any(diff(tt) <= 0)) {
    stop("tt must be finite, strictly increasing, and match U's rows.")
  }
  n_components <- ncol(U)
  attr(control, "state_bounds") <- validate_joint_gp_bounds(state_lower, state_upper, n_components, "state")
  n_parameters <- if (is.null(p0)) detect_n_params(f) else length(p0)
  if (n_parameters < 1L) stop("Supply p0 so the number of parameters is known.")
  if (is.null(p0)) p0 <- rep(1, n_parameters)
  if (!is.numeric(p0) || any(!is.finite(p0))) stop("p0 must be a finite numeric vector.")
  bounds <- validate_joint_gp_bounds(lower, upper, n_parameters, "parameter")
  lower <- bounds$lower
  upper <- bounds$upper
  p0 <- pmax(lower, pmin(upper, p0))
  if (!is.null(noise_sd)) {
    if (!is.numeric(noise_sd) || !length(noise_sd) %in% c(1L, n_components) || any(!is.finite(noise_sd) | noise_sd <= 0)) {
      stop("noise_sd must be NULL or a positive scalar/per-component vector.")
    }
    noise_sd <- rep_len(noise_sd, n_components)
  }
  if (control$verbose) message("Fitting ", n_components, " component GP(s).")
  fits <- lapply(seq_len(n_components), function(component) {
    fit_component_gp(
      tt, U[, component],
      if (is.null(noise_sd)) NULL else noise_sd[component], control
    )
  })
  if (control$weak_integration == "gp_gauss") {
    answer <- solve_with_gauss_quadrature(f, U, tt, p0, fits, control, lower, upper)
    return(answer)
  }
  model <- build_symbolic_ode_model(f, n_components, n_parameters, control$em_order)
  if (is.null(control$weak_radii) && control$weak_radius_method %in% c("sensitivity", "svd")) {
    answer <- solve_with_selected_weak_tests(U, tt, p0, fits, model, control, lower, upper)
    if (!answer$diagnostics$design_passed) {
      warning("No fully screened candidate space was available; using the reported diagnostic fallback.", call. = FALSE)
    }
    if (!answer$diagnostics$grid_passed) {
      warning_text <- paste0(
        "Observation-operator error exceeds grid_tol on the fixed working grid; ",
        "inspect interpolation diagnostics or explicitly choose a larger grid_min."
      )
      if (control$grid_action == "error") stop(warning_text) else warning(warning_text, call. = FALSE)
    }
    return(answer)
  }
  grid <- build_joint_gp_grid(tt, control)
  design <- build_weak_test_design(tt, fits, control)
  if (control$verbose) message("Preparing fixed grid with ", length(grid), " points.")
  state <- prepare_gp_state(fits, grid, control)
  H <- build_observation_matrix(grid, tt)
  interpolation_error <- vapply(fits, measure_gp_interpolation_error, numeric(1), tt = grid, obs = tt, H = H)
  weak <- build_grid_weak_operator(grid, design, control)
  initial_parameters <- initialize_ode_parameters(state$mean, p0, weak, model, state$Sigma, control, lower, upper)
  weak_value <- evaluate_weak_residual(state$mean, initial_parameters, weak, model)
  ode_covariance <- build_ode_metric(weak, n_components, control, weak_value$Ju, state$Sigma)
  metric <- prepare_ode_metric(ode_covariance, weak, U, tt, control)
  weights <- metric$weights
  residual_error <- if (control$lambda > 0) {
    check_initial_grid_refinement(state, initial_parameters, weak, model, fits, weights, control)
  } else {
    c(interior = 0, boundary = 0)
  }
  weak_problem <- list(weak = weak, Omega = ode_covariance, diagnostics = design$radius_selection)
  fit <- fit_observed_joint_gp(
    U, tt, fits, model, control, lower, upper,
    state, H, initial_parameters, weak_problem, residual_error
  )
  fit$diagnostics$initial_interpolation <- interpolation_error
  fit$diagnostics$grid_passed <- fit$diagnostics$grid_passed && max(interpolation_error) <= control$grid_tol
  if (!fit$diagnostics$grid_passed) {
    warning_text <- "Observation-operator error exceeds grid_tol. Inspect fit$diagnostics$final_grid."
    if (control$grid_action == "error") stop(warning_text) else warning(warning_text, call. = FALSE)
  }
  fit
}

fit_observed_joint_gp <- function(Y, tt, fits, model, control, lower, upper, state, H,
                                  pilot, weak_problem, initial_error) {
  weak <- weak_problem$weak
  n_components <- model$D
  n_parameters <- model$J
  noise <- sqrt(vapply(fits, `[[`, numeric(1), "noise2"))
  metric <- prepare_ode_metric(weak_problem$Omega, weak, Y, tt, control)
  weights <- metric$weights
  initial_coordinates <- c(pilot, pack_initial_gp_state(state))
  objective <- make_observed_gp_objective(Y, H, state, weak, model, weights, noise, control$lambda)
  coordinates <- state$coordinates
  coordinates$pscale <- rep(1, n_parameters)
  bounds <- build_optimizer_bounds(coordinates, lower, upper, state$bounds)
  result <- optimize_joint_gp(initial_coordinates, objective, control, bounds$lower, bounds$upper)
  point <- unpack_joint_gp_coordinates(result$par, coordinates)
  final <- objective(result$par, TRUE)
  stationarity <- diagnose_joint_gp_fit(point$U, point$p, Y, H, noise, weak, model, weights, control,
    state, if (state$include_gp_prior) seq_len(n_components) else integer(),
    lower = lower, upper = upper
  )
  accuracy <- check_joint_gp_accuracy(point$U, point$p, state, weak, model, fits, weights, H, control)
  initial <- objective(initial_coordinates, FALSE)
  mean_squared_weak_residual <- function(r) sum((weights$Wi %*% r)^2) / weights$interior$rank
  structure(list(
    phat = point$p, U_hat = point$U, U_obs_hat = H %*% point$U,
    weak_integration = control$weak_integration,
    tt = weak$tt, tt_obs = tt, Y = Y, gp = fits, noise_sd = noise,
    initial = list(p = pilot, U = state$mean, theta = initial_coordinates), lambda = control$lambda,
    ode_weighting = control$ode_weighting,
    ode_units = metric$scaling$units,
    include_gp_prior = state$include_gp_prior,
    objective = result$objective, contributions = final$contributions, converged = stationarity$stationary,
    optimizer = list(
      method = result$method, converged = result$converged, reason = result$reason,
      iterations = result$iterations, evaluations = result$evaluations
    ),
    convergence_reason = result$reason, iterations = result$iterations,
    diagnostics = list(
      rho_initial = mean_squared_weak_residual(initial$raw), rho_final = mean_squared_weak_residual(final$raw),
      ode_scaling = metric$scaling,
      state_coordinates = if (state$prior_whitened) "gp_whitened" else "physical_or_mixed",
      radius_selection = weak_problem$diagnostics, stationarity = stationarity,
      scaled_gradient = stationarity$scaled_gradient,
      design_passed = !identical(weak_problem$diagnostics$screening$passed, FALSE),
      interior_rank = weights$interior$rank, boundary_rank = weights$boundary$rank,
      gp_jitter = state$jitter, gp_converged = vapply(fits, function(g) g$convergence == 0, logical(1)),
      grid_passed = accuracy$passed,
      weak_grid_passed = max(initial_error) <= control$weak_grid_tol && accuracy$weak_passed,
      initial_error = initial_error, final_grid = accuracy
    ),
    problem = list(
      evaluate = objective, theta = result$par, state = state, H = H,
      weak = weak, model = model, Omega = metric$Omega, Omega_raw = weak_problem$Omega, weights = weights,
      control = control, lower = lower, upper = upper, state_bounds = state$bounds
    )
  ), class = "jointgp")
}

#' @export
print.jointgp <- function(x, ...) {
  cat("WENDyGP (", x$formulation,
    "; ", x$gp[[1]]$control$kernel, " covariance)\n",
    sep = ""
  )
  cat("Parameters:", format(x$phat, digits = 6), "\n")
  cat("Objective:", format(x$objective, digits = 6), " lambda:", x$lambda, "\n")
  cat("ODE weighting:", x$problem$control$ode_weighting, "\n")
  cat("ODE units:", x$ode_units, "\n")
  if (identical(x$weak_integration, "gp_gauss")) {
    cat("Weak integration: continuous GP / Gauss (fixed optimized grid)\n")
  }
  cat("GP prior:", if (isFALSE(x$include_gp_prior)) "disabled" else "enabled", "\n")
  cat("Converged:", x$converged, "\n")
  cat("Optimizer exit:", x$optimizer$reason, "\n")
  if (!x$converged) cat("Fit status:", x$convergence_reason, "\n")
  cat("Observation-operator accuracy passed:", x$diagnostics$grid_passed, "\n")
  invisible(x)
}
