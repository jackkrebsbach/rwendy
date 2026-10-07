# Latent-component estimation with a fixed weak operator and separate GP fit.
#
# Stage 1 fits one latent component in scaled physical coordinates, without a
# latent GP penalty. Stage 2 fits a GP to that curve. Stage 3 uses the fitted GP
# covariance to precondition the joint solve; latent_penalty optionally adds
# its quadratic penalty. Fitting the latent curve and its GP in separate stages
# avoids an unbounded joint objective over the state and GP amplitude.

identify_latent_component <- function(U) {
  observed <- colSums(is.na(U)) == 0L
  latent <- colSums(!is.na(U)) == 0L
  if (any(!observed & !latent)) {
    stop("Each component column must be fully observed or entirely NA.")
  }
  if (!any(observed)) stop("At least one component must be observed.")
  if (sum(latent) != 1L) stop("Exactly one latent component is supported.")
  if (any(!is.finite(U[, observed, drop = FALSE]))) stop("Observed components must be finite.")
  list(observed = which(observed), latent = which(latent))
}

# Initialize the latent component from the observed row mean.
initialize_latent_state <- function(U, tt, grid, observed, latent, latent_initial = NULL) {
  initial_state <- matrix(0, length(grid), ncol(U))
  for (component in observed) {
    initial_state[, component] <- stats::approx(tt, U[, component], xout = grid, rule = 2)$y
  }
  curve <- rowMeans(initial_state[, observed, drop = FALSE])
  observed_rms <- sqrt(colMeans(initial_state[, observed, drop = FALSE]^2))
  if (sqrt(mean(curve^2)) <= sqrt(.Machine$double.eps) * max(observed_rms)) {
    curve <- initial_state[, observed[which.max(observed_rms)]]
  }
  if (!is.null(latent_initial)) {
    curve <- stats::approx(tt, latent_initial, xout = grid, rule = 2)$y
  }
  initial_state[, latent] <- curve
  initial_state
}

# GP posterior mean at arbitrary times.
predict_gp_mean <- function(fit, tt) {
  as.vector(fit$mean + gp_cross_covariance(fit, tt, fit$tt) %*% fit$alpha)
}

calculate_latent_reference_scale <- function(initial_state, latent) {
  lapply(latent, function(component) {
    scale <- max(stats::sd(initial_state[, component]), .1 * sqrt(mean(initial_state[, component]^2)), 1e-6)
    list(scale = scale, mean = mean(initial_state[, component]) / scale)
  })
}

# Stationary Matern 5/2 at a fixed fraction of the normalized domain. Radius
# bounds are grid/domain based, never chosen using latent truth.
build_latent_interpolation_kernel <- function(problem, radius = 1 / 16) {
  span <- diff(range(problem$tt))
  bounds <- c(2 * mean(diff(problem$grid)) / span, 2)
  radius_coefficient <- stats::qlogis(pmin(.99, pmax(.01, log(radius / bounds[1]) / diff(log(bounds)))))
  lapply(seq_along(problem$latent), function(k) {
    c(problem$reference[[k]], list(
      origin = min(problem$tt), span = span,
      bounds = bounds, radius_coefficient = radius_coefficient
    ))
  })
}

# Only observed components have a prior before stage 2 fits the latent GP.
prepare_observed_gp_prior <- function(problem) {
  n_grid <- length(problem$grid)
  n_components <- problem$n_components
  prior_mean <- matrix(0, n_grid, n_components)
  prior_factors <- vector("list", n_components)
  for (k in seq_along(problem$observed_fits)) {
    component <- problem$observed[k]
    fit <- problem$observed_fits[[k]]
    prior_mean[, component] <- fit$mean
    covariance <- gp_cross_covariance(fit, problem$grid, problem$grid)
    prior_factors[[component]] <- t(factor_gp_covariance(covariance, problem$control$gp_jitter)$R)
  }
  list(mu = prior_mean, L = prior_factors)
}

# Fixed posterior coordinates precondition the observed data+GP block. The
# objective still contains the original observations and the original prior;
# only the coordinates the optimizer moves in change.
prepare_latent_coordinates <- function(problem) {
  n_grid <- length(problem$grid)
  prior_state <- prepare_observed_gp_prior(problem)
  observed_coordinates <- prior_state
  for (k in seq_along(problem$observed)) {
    component <- problem$observed[k]
    L <- observed_coordinates$L[[component]]
    A <- problem$H %*% L / problem$noise[k]
    R <- chol(diag(n_grid) + crossprod(A))
    R_inv <- backsolve(R, diag(n_grid))
    observed_coordinates$L[[component]] <- t(chol(tcrossprod(L %*% R_inv)))
    observed_coordinates$mu[, component] <- prior_state$mu[, component] + L %*% solve_cholesky_system(
      R,
      crossprod(A, (problem$Y[, k] - problem$H %*% prior_state$mu[, component]) / problem$noise[k])
    )
  }
  list(
    n_grid = n_grid,
    latent_component = problem$latent[1],
    observed_coordinates = observed_coordinates,
    prior_state = prior_state,
    latent_scale = problem$reference[[1]]$scale,
    parameter_scale = pmax(abs(problem$p0), .1)
  )
}

# Order the component fits and prior factors for quadrature. The initial latent
# kernel only defines interpolation; it contributes no prior penalty.
prepare_latent_weak_extension <- function(problem) {
  n_components <- problem$n_components
  kernels <- build_latent_interpolation_kernel(problem)
  fits <- vector("list", n_components)
  for (k in seq_along(problem$observed)) {
    fits[[problem$observed[k]]] <- problem$observed_fits[[k]]
  }
  for (k in seq_along(problem$latent)) {
    component <- problem$latent[k]
    kernel <- kernels[[k]]
    fits[[component]] <- list(
      origin = kernel$origin, span = kernel$span,
      tau2 = kernel$scale^2, radius_coef = kernel$radius_coefficient,
      radius_bounds = kernel$bounds, control = problem$control,
      mean = kernel$scale * kernel$mean, noise2 = problem$noise[1]^2
    )
  }
  prior_factors <- lapply(fits, function(f) {
    covariance <- gp_cross_covariance(f, problem$grid, problem$grid)
    t(factor_gp_covariance(covariance, problem$control$gp_jitter)$R)
  })
  list(fits = fits, state = list(L = prior_factors))
}

build_latent_weak_problem <- function(problem) {
  control <- problem$control
  n_components <- problem$n_components
  extension <- prepare_latent_weak_extension(problem)
  weak <- build_gauss_weak_operator(problem$grid, problem$design, extension$fits, extension$state, control, problem$maps)
  weak$parameter_scale <- if (is.null(control$weak_design_scale)) pmax(abs(problem$p0), .1) else control$weak_design_scale
  if (length(weak$parameter_scale) != problem$n_parameters) stop("weak_design_scale must have one value per parameter.")
  weak$state_scale <- problem$units
  metric <- build_ode_metric(weak, n_components, control)
  weights <- build_ode_weights(metric, weak, n_components, control)
  problem$weak <- weak
  problem$weights <- scale_weak_weights(weights, problem$scaling$precision, weak$K)
  problem$extension_fits <- extension$fits
  problem$extension_state <- extension$state
  problem
}

# Stage 1 uses observed posterior coordinates and scaled physical latent values.
make_initial_latent_objective <- function(problem, coordinate_setup) {
  coordinates <- coordinate_setup$observed_coordinates
  coordinates$mu[, coordinate_setup$latent_component] <- 0
  coordinates$L[[coordinate_setup$latent_component]] <- coordinate_setup$latent_scale
  coordinates$order <- c(problem$observed, coordinate_setup$latent_component)
  coordinates$pscale <- coordinate_setup$parameter_scale
  make_latent_gp_objective(problem, coordinate_setup, coordinates)
}

# Both latent stages use exactly the same residual engine as observed fits.
make_latent_gp_objective <- function(problem, coordinate_setup, coordinates, prior = NULL, penalty = FALSE) {
  coordinates <- apply_state_bounds_to_coordinates(coordinates, attr(problem$control, "state_bounds"), problem$units)
  state <- coordinate_setup$prior_state
  components <- problem$observed
  if (penalty) {
    state$L[[coordinate_setup$latent_component]] <- prior$L
    state$mu[, coordinate_setup$latent_component] <- prior$mean
    components <- c(components, coordinate_setup$latent_component)
  }
  evaluate <- make_joint_gp_objective(
    problem$Y, problem$H, problem$weak, problem$model, problem$weights,
    problem$noise, problem$control$lambda, coordinates, state, components, problem$observed
  )
  objective <- function(par, jacobian = TRUE) {
    value <- evaluate(par, jacobian)
    contributions <- value$contributions
    value$contributions <- c(
      data = unname(contributions["data"]),
      observed_gp = sum(value$prior_contributions[as.character(problem$observed)]),
      ode = sum(contributions[c("interior", "boundary")])
    )
    if (penalty) {
      value$contributions <- c(value$contributions,
        latent_gp = unname(value$prior_contributions[as.character(coordinate_setup$latent_component)])
      )
    }
    value
  }
  attr(objective, "coordinates") <- coordinates
  attr(objective, "bounds") <- build_optimizer_bounds(
    coordinates, problem$lower, problem$upper,
    attr(problem$control, "state_bounds")
  )
  objective
}

fit_initial_latent_state <- function(problem, coordinate_setup, multiplier, maxit = 5000L,
                                     gradient_tol = 1e-4) {
  n_parameters <- length(problem$p0)
  n_grid <- coordinate_setup$n_grid
  n_observed_coordinates <- n_grid * length(problem$observed)
  U <- problem$pilot
  p <- problem$p0 * multiplier
  objective <- make_initial_latent_objective(problem, coordinate_setup)
  bounds <- attr(objective, "bounds")
  par <- pack_joint_gp_coordinates(p, U, attr(objective, "coordinates"))
  par <- pmax(bounds$lower, pmin(bounds$upper, par))
  # Warm phase on the parameters and the latent only, holding the observed
  # blocks at their posterior-preconditioned start.
  warm <- optimize_joint_gp(par, objective, problem$control,
    maxit = min(75L, maxit),
    active = c(seq_len(n_parameters), n_parameters + n_observed_coordinates + seq_len(n_grid)),
    lower = bounds$lower, upper = bounds$upper
  )
  par <- warm$par
  fit <- optimize_joint_gp(par, objective, problem$control, bounds$lower, bounds$upper, maxit = maxit)
  value <- objective(fit$par)
  stationarity <- diagnose_latent_fit(problem, coordinate_setup, value$p, value$U, gradient_tol)
  list(
    value = value$value, p = value$p, U = value$U,
    contributions = value$contributions, multiplier = multiplier, theta = fit$par,
    converged = stationarity$stationary, stationarity = stationarity,
    optimizer = fit[setdiff(names(fit), "par")], warm = warm[setdiff(names(warm), "par")],
    projected_gradient = max(abs(value$gradient))
  )
}

diagnose_latent_fit <- function(problem, coordinate_setup, p, U, tolerance, prior = NULL, penalty = FALSE, rank = FALSE) {
  state <- coordinate_setup$prior_state
  components <- problem$observed
  if (penalty) {
    state$L[[coordinate_setup$latent_component]] <- prior$L
    state$mu[, coordinate_setup$latent_component] <- prior$mean
    components <- c(components, coordinate_setup$latent_component)
  }
  diagnose_joint_gp_fit(U, p, problem$Y, problem$H, problem$noise, problem$weak, problem$model,
    problem$weights, problem$control, state, components,
    observed = problem$observed,
    lower = problem$lower, upper = problem$upper, tolerance = tolerance, rank = rank
  )
}

# Stage 2 fits the latent GP to the initial curve. Its nugget absorbs roughness
# in that curve and is excluded from the prior covariance.
fit_latent_gp_prior <- function(problem, u_latent) {
  fit <- fit_component_gp(problem$grid, u_latent, NULL, problem$control)
  covariance <- gp_cross_covariance(fit, problem$grid, problem$grid)
  covariance <- (covariance + t(covariance)) / 2
  factor <- factor_gp_covariance(covariance, problem$control$gp_jitter)
  list(
    fit = fit, mean = fit$mean, tau = sqrt(fit$tau2), L = t(factor$R),
    smooth = predict_component_gp(fit, problem$grid)$mean,
    nugget_sd = sqrt(fit$noise2), jitter = factor$added
  )
}

# Stage 3 keeps the fitted latent GP fixed during the joint solve.
make_final_latent_objective <- function(problem, coordinate_setup, prior, penalty = FALSE) {
  coordinates <- coordinate_setup$observed_coordinates
  coordinates$L[[coordinate_setup$latent_component]] <- prior$L
  coordinates$mu[, coordinate_setup$latent_component] <- prior$mean
  coordinates$order <- seq_len(problem$n_components)
  coordinates$pscale <- coordinate_setup$parameter_scale
  make_latent_gp_objective(problem, coordinate_setup, coordinates, prior, penalty)
}

fit_final_latent_state <- function(problem, coordinate_setup, prior, p0, U0, maxit = 5000L,
                                   penalty = FALSE, gradient_tol = 1e-4) {
  objective <- make_final_latent_objective(problem, coordinate_setup, prior, penalty)
  bounds <- attr(objective, "bounds")
  par <- pack_joint_gp_coordinates(p0, U0, attr(objective, "coordinates"))
  fit <- optimize_joint_gp(par, objective, problem$control, bounds$lower, bounds$upper, maxit = maxit)
  value <- objective(fit$par)
  stationarity <- diagnose_latent_fit(problem, coordinate_setup, value$p, value$U, gradient_tol, prior, penalty)
  list(
    value = value$value, p = value$p, U = value$U, contributions = value$contributions,
    latent_z = sum(forwardsolve(prior$L, value$U[, coordinate_setup$latent_component] - prior$mean)^2) / 2,
    converged = stationarity$stationary, stationarity = stationarity,
    theta = fit$par,
    optimizer = fit[setdiff(names(fit), "par")],
    projected_gradient = max(abs(value$gradient))
  )
}

# Full weak test basis and fixed component units.
prepare_latent_problem <- function(f, U, tt, p0, lower, upper, noise_sd, control,
                                   latent_initial = NULL) {
  components <- identify_latent_component(U)
  observed <- components$observed
  latent <- components$latent
  n_components <- ncol(U)
  n_parameters <- length(p0)
  Y <- U[, observed, drop = FALSE]
  observed_noise_sd <- if (is.null(noise_sd)) NULL else rep_len(noise_sd, length(observed))
  observed_fits <- lapply(seq_along(observed), function(k) {
    fit_component_gp(
      tt, Y[, k],
      if (is.null(observed_noise_sd)) NULL else observed_noise_sd[k], control
    )
  })
  grid <- build_joint_gp_grid(tt, control)
  model <- build_symbolic_ode_model(f, n_components, n_parameters, 0L)
  scale_pilot <- initialize_latent_state(U, tt, grid, observed, latent)
  initial_state <- initialize_latent_state(U, tt, grid, observed, latent, latent_initial)
  bounds <- attr(control, "state_bounds")
  if (!is.null(bounds)) {
    # Smooth bounded observed starts before interpolation: high-order extension
    # of raw noise can leave a model's domain even when its grid values do not.
    for (k in seq_along(observed)) {
      if (bounds$bounded[observed[k]]) {
        initial_state[, observed[k]] <- predict_gp_mean(observed_fits[[k]], grid)
      }
    }
    initial_state <- clip_state_to_bounds(initial_state, bounds)
  }
  design <- if (is.null(control$weak_radii)) {
    build_weak_test_pool(grid, control)
  } else {
    build_weak_test_design(tt, observed_fits, control)
  }
  maps <- list(
    interior = orthonormalize_weak_tests(evaluate_weak_test_rows(grid, design$interior, 0L, control$bump_eta), control$basis_tol),
    boundary = orthonormalize_weak_tests(evaluate_weak_test_rows(grid, design$boundary, 0L, control$bump_eta), control$basis_tol)
  )
  # Freeze the default scale policy independently of a supplied initial curve.
  # This remains a heuristic for unknown latent units; explicit component
  # scales are supported and reported for comparisons across physical units.
  units <- sqrt(colMeans(scale_pilot^2))
  units[observed] <- sqrt(colMeans(Y^2))
  if (any(!is.finite(units) | units <= 0)) stop("Unidentified pilot component scale.")
  scaling <- calculate_ode_scaling(matrix(rep(units, each = length(tt)), length(tt)), tt, control)
  scaling$scale_source <- if (is.null(control$ode_component_scale)) {
    "observed_rms_and_latent_initializer_heuristic"
  } else {
    "supplied"
  }
  units <- scaling$rms
  problem <- list(
    Y = Y, tt = tt, grid = grid, H = build_observation_matrix(grid, tt), observed_fits = observed_fits,
    noise = sqrt(vapply(observed_fits, `[[`, numeric(1), "noise2")), model = model,
    design = design, maps = maps, units = units, scaling = scaling, pilot = initial_state,
    n_components = n_components, n_parameters = n_parameters,
    observed = observed, latent = latent, control = control, p0 = p0, lower = lower,
    upper = upper
  )
  problem$reference <- calculate_latent_reference_scale(initial_state, latent)
  problem
}

#' Compatibility wrapper for latent WENDyGP estimation
#'
#' @description Calls [solveWendyGP()] with formulation="latent". New code
#'   should use that unified entry point.
#'
#' @param f Function f(u,p,t), compatible with WENDy's symbolic machinery.
#' @param U Numeric observation matrix with exactly one entirely NA column.
#'   Every observed column must be fully finite, with at least six observations.
#' @param tt Strictly increasing finite observation times.
#' @param p0 Initial parameters, or NULL to infer their count and start at ones.
#' @param noise_sd Positive scalar or per-observed-component noise SD, or NULL.
#' @param control Overrides from [wendygp_control()]. This formulation requires
#'   matern52, gp_gauss, observed GP priors, and test_gram or identity weighting.
#'   It uses the full numerical test span. weak_radii and weak_design_radii
#'   control explicit tests and the multiscale pool respectively. The observed
#'   extension is Lagrange and the latent extension is GP. ode_units and
#'   ode_component_scale control the ODE metric. Supply component scales when
#'   the latent units cannot be inferred from the observed components; the
#'   default latent scale is an initializer-based heuristic, not unit invariant.
#' @param starts Finite parameter multipliers. Stationary completed runs are
#'   preferred, followed by finite feasible stage-3 estimates, then stage-1
#'   fallbacks. Runs in the preferred group are ranked by objective value;
#'   all runs are retained. With latent_penalty=TRUE
#'   only one start is supported, so different fitted priors are not compared.
#' @param latent_penalty Add a frozen latent GP quadratic in stage 3. FALSE uses
#'   its fitted covariance only to change optimization coordinates for unbounded
#'   states and starts from the stage-1 estimate. Bounded states use component
#'   scaling. TRUE starts
#'   from the smoothed stage-1 curve for the new penalized objective.
#' @param maxit Iteration cap per optimization phase, default 5000. An explicit
#'   argument takes precedence over control$maxit; otherwise that override is used.
#' @param gradient_tol Maximum Jacobian-column-scaled projected gradient in
#'   physical coordinates, for stages 1 and 3. Defaults to 1e-4; an explicit
#'   argument takes precedence over control$gtol. No objective normalization is used.
#' @param latent_initial Optional finite latent starting curve at observation
#'   times. It does not change the frozen component scales. By default the row
#'   mean of observations is used, falling back to the observed column of
#'   largest RMS if the row mean cancels numerically.
#' @inheritParams solveWendyGP
#' @return A jointgp object (also inheriting wendygp_latent) with parameter/state estimates, all start
#'   records, optimizer exits, physical stationarity, local Jacobian
#'   rank, and separate grid-resolution and quadrature diagnostics. A finite
#'   feasible stage-3 estimate is returned even without stationarity, with a
#'   warning and converged=FALSE. Stage 1 is returned only if no usable stage-3
#'   estimate exists. Numerical errors stop the solve.
#'   latent_prior contains the separately fitted GP, including its nugget.
#'   Local stationarity and rank do not establish global identifiability.
#' @export
solveWendyGPLatent <- function(f, U, tt, p0 = NULL, noise_sd = NULL, control = NULL,
                               starts = 1, latent_penalty = FALSE,
                               maxit = 5000L, gradient_tol = 1e-4,
                               latent_initial = NULL,
                               parameter_lower = NULL, parameter_upper = NULL,
                               state_lower = NULL, state_upper = NULL, lower = NULL, upper = NULL) {
  # Compatibility only: validation, preparation and fitting live in solveWendyGP.
  args <- list(
    f = f, U = U, tt = tt, p0 = p0, noise_sd = noise_sd, control = control,
    formulation = "latent", starts = starts,
    latent_penalty = latent_penalty, latent_initial = latent_initial,
    parameter_lower = parameter_lower, parameter_upper = parameter_upper,
    lower = lower, upper = upper, state_lower = state_lower, state_upper = state_upper
  )
  if (!missing(maxit)) args$maxit <- maxit
  if (!missing(gradient_tol)) args$gradient_tol <- gradient_tol
  answer <- do.call(solveWendyGP, args)
  answer$call <- match.call()
  answer
}

solve_latent_joint_gp <- function(f, U, tt, p0 = NULL, noise_sd = NULL, control = NULL,
                                  starts = 1,
                                  latent_penalty = FALSE, maxit = 5000L,
                                  gradient_tol = 1e-4,
                                  latent_initial = NULL,
                                  lower = NULL, upper = NULL, state_lower = NULL, state_upper = NULL) {
  supplied <- if (is.null(control)) list() else control
  # Validate unsupported choices before the observed control resolver can
  # coerce them to another integration scheme.
  if (!is.null(supplied$kernel) && !identical(supplied$kernel, "matern52")) {
    stop("Latent estimation currently supports kernel='matern52' only.")
  }
  if (isFALSE(supplied$include_gp_prior)) stop("Latent estimation requires observed GP priors.")
  if (identical(supplied$ode_weighting, "gp_delta")) {
    stop("Latent estimation supports ode_weighting='test_gram' or 'identity', not 'gp_delta'.")
  }
  if (!is.null(supplied$weak_radius_method) && supplied$weak_radius_method != "svd") {
    stop("Latent estimation requires weak_radius_method='svd'; use weak_radii for explicit tests.")
  }
  control <- do.call(wendygp_control, supplied)
  if (control$weak_integration != "gp_gauss") stop("Latent estimation requires weak_integration='gp_gauss'.")
  if (!is.null(control$weak_design_budget)) stop("Latent estimation uses the full test span; weak_design_budget is unsupported.")
  if (!is.matrix(U) || !is.numeric(U) || nrow(U) < 6L || ncol(U) < 2L) {
    stop("U must be a numeric matrix with at least six rows and two components.")
  }
  if (!is.numeric(tt) || length(tt) != nrow(U) || any(!is.finite(tt)) || any(diff(tt) <= 0)) {
    stop("tt must be strictly increasing, finite, and match nrow(U).")
  }
  components <- identify_latent_component(U)
  n_components <- ncol(U)
  if (!is.numeric(starts) || !length(starts) || any(!is.finite(starts))) {
    stop("starts must contain finite parameter multipliers.")
  }
  if (!is.logical(latent_penalty) || length(latent_penalty) != 1L || is.na(latent_penalty)) {
    stop("latent_penalty must be TRUE or FALSE.")
  }
  if (latent_penalty && length(starts) > 1L) {
    stop("latent_penalty=TRUE requires one start; different stage-2 priors cannot be ranked as one fixed objective.")
  }
  if (!is.null(noise_sd) && (!is.numeric(noise_sd) || !length(noise_sd) %in% c(1L, length(components$observed)) ||
    any(!is.finite(noise_sd)) || any(noise_sd <= 0))) {
    stop("noise_sd must be positive, scalar or per observed component.")
  }
  if (!is.null(latent_initial) && (!is.numeric(latent_initial) ||
    length(latent_initial) != nrow(U) || any(!is.finite(latent_initial)))) {
    stop("latent_initial must be finite with one value per observation time.")
  }
  n_parameters <- if (is.null(p0)) detect_n_params(f) else length(p0)
  if (n_parameters < 1L) stop("Supply p0 so the number of parameters is known.")
  if (is.null(p0)) p0 <- rep(1, n_parameters)
  if (!is.numeric(p0) || any(!is.finite(p0))) stop("p0 must be a finite numeric vector.")
  bounds <- validate_joint_gp_bounds(lower, upper, n_parameters, "parameter")
  lower <- bounds$lower
  upper <- bounds$upper
  p0 <- pmax(lower, pmin(upper, p0))
  attr(control, "state_bounds") <- validate_joint_gp_bounds(state_lower, state_upper, n_components, "state")
  expected_extension <- ifelse(seq_len(n_components) %in% components$observed, "lagrange", "gp")
  if ("weak_extension" %in% names(supplied) &&
    !identical(supplied$weak_extension, "lagrange") &&
    !identical(supplied$weak_extension, expected_extension)) {
    stop("The latent formulation uses Lagrange for observed components and GP for the latent component.")
  }
  control$weak_extension <- expected_extension
  control$em_order <- 0L
  control$weak_design_info <- 1
  control$maxit <- maxit
  control$gtol <- gradient_tol
  initial_problem <- prepare_latent_problem(f, U, tt, p0, lower, upper, noise_sd, control, latent_initial)
  problem <- build_latent_weak_problem(initial_problem)
  coordinate_setup <- prepare_latent_coordinates(problem)
  accuracy <- function(U, p) {
    check_gauss_quadrature_accuracy(U, p, problem$extension_state, problem$weak,
      problem$model, problem$extension_fits, problem$weights, problem$H, control,
      observed = problem$observed, observation_times = tt
    )
  }
  initial_accuracy <- accuracy(problem$pilot, p0)
  runs <- lapply(starts, function(multiplier) {
    run <- list(multiplier = multiplier, stage1 = NULL, prior = NULL, stage3 = NULL)
    run$stage1 <- fit_initial_latent_state(problem, coordinate_setup, multiplier, maxit, gradient_tol)
    run$stage1$accuracy <- accuracy(run$stage1$U, run$stage1$p)
    if (!run$stage1$converged) {
      run$status <- "not_stationary"
      run$message <- "Stage 1 did not reach physical stationarity."
      return(run)
    }
    run$prior <- fit_latent_gp_prior(problem, run$stage1$U[, coordinate_setup$latent_component])
    if (run$prior$fit$convergence != 0L) {
      run$status <- "gp_not_converged"
      run$message <- "Stage-2 GP fitting did not converge."
      return(run)
    }
    U0 <- run$stage1$U
    if (latent_penalty) U0[, coordinate_setup$latent_component] <- run$prior$smooth
    run$stage3 <- fit_final_latent_state(problem, coordinate_setup, run$prior, run$stage1$p, U0, maxit,
      penalty = latent_penalty, gradient_tol = gradient_tol
    )
    if (!run$stage3$converged) {
      run$status <- "not_stationary"
      run$message <- "Stage 3 did not reach physical stationarity."
      return(run)
    }
    run$status <- "stationary"
    run
  })
  usable <- function(x) {
    !is.null(x) && is.finite(x$value) &&
      all(is.finite(c(x$p, x$U))) && isTRUE(x$stationarity$feasible)
  }
  valid <- which(vapply(runs, function(r) identical(r$status, "stationary"), logical(1)))
  selected_stage <- vapply(runs, function(r) {
    if (usable(r$stage3)) 3L else if (usable(r$stage1)) 1L else 0L
  }, integer(1))
  # Return the final optimization's estimate even when it is nonstationary or
  # has a higher objective than stage 1. Keep stage 1 as a numerical fallback.
  stage3_available <- which(selected_stage == 3L)
  eligible <- if (length(valid)) {
    valid
  } else if (length(stage3_available)) {
    stage3_available
  } else {
    which(selected_stage == 1L)
  }
  if (!length(eligible)) stop("No finite, feasible latent estimate was produced.")
  values <- vapply(eligible, function(k) runs[[k]][[paste0("stage", selected_stage[k])]]$value, numeric(1))
  chosen <- eligible[which.min(values)]
  best <- runs[[chosen]]
  final_stage <- selected_stage[chosen]
  final <- best[[paste0("stage", final_stage)]]
  applied_penalty <- latent_penalty && final_stage == 3L
  stationarity <- diagnose_latent_fit(problem, coordinate_setup, final$p, final$U,
    gradient_tol, best$prior, applied_penalty,
    rank = TRUE
  )
  final_accuracy <- accuracy(final$U, final$p)
  pipeline_converged <- identical(best$status, "stationary") && stationarity$stationary
  reason <- if (pipeline_converged) {
    final$optimizer$reason
  } else {
    paste0(best$message, " Returning the stage-", final_stage, " estimate.")
  }
  if (!pipeline_converged) warning(reason, " See fit$runs; converged=FALSE.", call. = FALSE)
  grid_passed <- initial_accuracy$passed && best$stage1$accuracy$passed && final_accuracy$passed
  weak_passed <- final_accuracy$quadrature_passed
  if (!grid_passed) {
    warning_text <- "Latent fit observation-operator error exceeds grid_tol; inspect initial_quadrature and final_grid."
    if (control$grid_action == "error") stop(warning_text) else warning(warning_text, call. = FALSE)
  }
  if (!weak_passed) warning(format_quadrature_warning(final_accuracy, control), call. = FALSE)
  evaluate <- if (final_stage == 3L) {
    make_final_latent_objective(problem, coordinate_setup, best$prior, applied_penalty)
  } else {
    make_initial_latent_objective(problem, coordinate_setup)
  }
  structure(list(
    phat = final$p, U_hat = final$U,
    U_obs_hat = problem$H %*% final$U[, problem$observed, drop = FALSE],
    U_stage1 = best$stage1$U, U_smooth = best$prior$smooth,
    observed = problem$observed, latent = problem$latent, tt = problem$grid, tt_obs = tt,
    Y = problem$Y, noise_sd = problem$noise, gp = problem$observed_fits, latent_prior = best$prior,
    lambda = control$lambda, ode_weighting = control$ode_weighting, ode_units = control$ode_units,
    objective = final$value, converged = pipeline_converged,
    optimizer = final$optimizer, contributions = final$contributions,
    runs = runs, multiplier = best$multiplier, selected_start = chosen,
    diagnostics = list(
      final_stage = final_stage, pipeline_converged = pipeline_converged,
      stage_status = best$status, modes = nrow(problem$maps$interior),
      weak_rows = problem$weak$K * n_components, state_variables = length(final$U),
      stage1_variables = n_parameters + length(best$stage1$U),
      weak_extension = control$weak_extension, ode_scaling = problem$scaling,
      stationarity = stationarity, scaled_gradient = stationarity$scaled_gradient, gradient_tol = gradient_tol,
      grid_passed = grid_passed, weak_grid_passed = weak_passed,
      quadrature_passed = final_accuracy$quadrature_passed,
      initial_quadrature = initial_accuracy, stage1_grid = best$stage1$accuracy, final_grid = final_accuracy,
      gp_converged = vapply(problem$observed_fits, function(g) g$convergence == 0, logical(1)),
      stage1_value = best$stage1$value, stage1_converged = best$stage1$converged,
      latent_nugget_sd = best$prior$nugget_sd, latent_tau = best$prior$tau,
      latent_penalty = applied_penalty, latent_penalty_requested = latent_penalty, latent_z = final$latent_z
    ),
    problem = list(
      evaluate = evaluate, theta = final$theta,
      coordinates = attr(evaluate, "coordinates"), H = problem$H, weak = problem$weak, model = problem$model,
      weights = problem$weights, control = control, lower = lower, upper = upper,
      state_bounds = attr(control, "state_bounds")
    ),
    include_gp_prior = TRUE, weak_integration = "gp_gauss",
    convergence_reason = reason, iterations = final$optimizer$iterations,
    control = control, call = match.call()
  ), class = c("jointgp", "wendygp_latent"))
}
