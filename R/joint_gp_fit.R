# Joint fitting: bounds and optimizer coordinates, the stacked objective,
# projected Levenberg-Marquardt, then the observed and latent formulations.

# ---- Bounds and optimizer coordinates --------------------------------------
# Coordinates U = mu + L z. Bounded columns use a diagonal L, so box bounds on
# z are box bounds on the state.

validate_joint_gp_bounds <- function(lower, upper, size, label) {
  expand <- function(x, default, side) {
    if (is.null(x)) {
      return(rep(default, size))
    }
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x) %in% c(1L, size) || anyNA(x)) {
      stop(label, " ", side, " must be a numeric scalar or a vector of length ", size, " without NA.")
    }
    rep_len(x, size)
  }
  lower_bounds <- expand(lower, -Inf, "lower bound")
  upper_bounds <- expand(upper, Inf, "upper bound")
  if (any(lower_bounds >= upper_bounds)) stop("Each ", label, " lower bound must be smaller than its upper bound.")
  list(lower = lower_bounds, upper = upper_bounds, bounded = is.finite(lower_bounds) | is.finite(upper_bounds))
}

apply_state_bounds_to_coordinates <- function(coordinates, bounds, scales = rep(1, ncol(coordinates$mu))) {
  n_grid <- nrow(coordinates$mu)
  for (component in which(bounds$bounded)) {
    coordinates$mu[, component] <- 0
    coordinates$L[[component]] <- diag(scales[component], n_grid)
  }
  coordinates
}

build_optimizer_bounds <- function(coordinates, lower, upper, state_bounds) {
  n_components <- ncol(coordinates$mu)
  n_grid <- nrow(coordinates$mu)
  state_lower <- rep(-Inf, n_grid * n_components)
  state_upper <- rep(Inf, n_grid * n_components)
  for (component in which(state_bounds$bounded)) {
    scale <- diag(coordinates$L[[component]])
    indices <- (match(component, coordinates$order) - 1L) * n_grid + seq_len(n_grid)
    state_lower[indices] <- (state_bounds$lower[component] - coordinates$mu[, component]) / scale
    state_upper[indices] <- (state_bounds$upper[component] - coordinates$mu[, component]) / scale
  }
  list(lower = c(lower / coordinates$pscale, state_lower), upper = c(upper / coordinates$pscale, state_upper))
}

pack_joint_gp_coordinates <- function(p, U, coordinates) {
  c(p / coordinates$pscale, unlist(lapply(coordinates$order, function(component) {
    L <- coordinates$L[[component]]
    centered <- U[, component] - coordinates$mu[, component]
    # BLAS solves multiply by the reciprocal; divide diagonal maps exactly.
    if (all(L[lower.tri(L)] == 0)) centered / diag(L) else forwardsolve(L, centered)
  })))
}

unpack_joint_gp_coordinates <- function(theta, coordinates) {
  n_parameters <- length(coordinates$pscale)
  n_grid <- nrow(coordinates$mu)
  U <- coordinates$mu
  for (k in seq_along(coordinates$order)) {
    component <- coordinates$order[k]
    state_coordinates <- theta[n_parameters + (k - 1L) * n_grid + seq_len(n_grid)]
    U[, component] <- U[, component] + as.vector(coordinates$L[[component]] %*% state_coordinates)
  }
  list(p = theta[seq_len(n_parameters)] * coordinates$pscale, U = U)
}

clip_state_to_bounds <- function(U, bounds) {
  for (component in which(bounds$bounded)) {
    U[, component] <- pmax(bounds$lower[component], pmin(bounds$upper[component], U[, component]))
  }
  U
}

# ---- Objective and optimizer -----------------------------------------------
# Stacked residual [data; GP prior; sqrt(lambda) W R] and its Jacobian in
# optimizer coordinates. Only prior_components get GP prior rows.
make_joint_gp_objective <- function(Y, H, weak, model, weights, noise, lambda,
                                    coordinates, prior_state = NULL,
                                    prior_components = integer(),
                                    observed = seq_len(model$D)) {
  n_grid <- nrow(coordinates$mu)
  n_components <- model$D
  n_parameters <- model$J
  n_observations <- nrow(Y)
  state_columns <- lapply(seq_len(n_components), function(component) {
    n_parameters + (match(component, coordinates$order) - 1L) * n_grid + seq_len(n_grid)
  })
  transform_state_jacobian <- function(block, component) block %*% coordinates$L[[component]]

  # Data and prior rows are linear, so their Jacobian is fixed.
  data_rows <- n_observations * length(observed)
  prior_rows <- n_grid * length(prior_components)
  linear_jacobian <- matrix(0, data_rows + prior_rows, n_parameters + n_grid * n_components)
  for (k in seq_along(observed)) {
    component <- observed[k]
    rows <- (k - 1L) * n_observations + seq_len(n_observations)
    linear_jacobian[rows, state_columns[[component]]] <- -transform_state_jacobian(H, component) / noise[k]
  }
  for (k in seq_along(prior_components)) {
    component <- prior_components[k]
    rows <- data_rows + (k - 1L) * n_grid + seq_len(n_grid)
    linear_jacobian[rows, state_columns[[component]]] <-
      forwardsolve(prior_state$L[[component]], transform_state_jacobian(diag(n_grid), component))
  }

  function(theta, jacobian = TRUE) {
    point <- unpack_joint_gp_coordinates(theta, coordinates)
    U <- point$U
    p <- point$p
    data_residual <- sweep(Y - H %*% U[, observed, drop = FALSE], 2, noise, "/")
    prior_residual <- matrix(0, n_grid, length(prior_components))
    for (k in seq_along(prior_components)) {
      component <- prior_components[k]
      prior_residual[, k] <- forwardsolve(prior_state$L[[component]], U[, component] - prior_state$mu[, component])
    }
    weak_value <- evaluate_weak_residual(U, p, weak, model, jacobian)
    ode_residual <- as.vector(sqrt(lambda) * weights$W %*% weak_value$r)
    residual <- c(data_residual, prior_residual, ode_residual)
    interior <- seq_len(nrow(weights$Wi))
    boundary <- length(interior) + seq_len(nrow(weights$Wb))
    contributions <- c(
      data = sum(data_residual^2) / 2,
      gp = sum(prior_residual^2) / 2,
      interior = sum(ode_residual[interior]^2) / 2,
      boundary = sum(ode_residual[boundary]^2) / 2
    )

    residual_jacobian <- NULL
    gradient <- NULL
    if (jacobian) {
      state_jacobian <- matrix(0, nrow(weak_value$Ju), n_grid * n_components)
      for (component in seq_len(n_components)) {
        columns <- (component - 1L) * n_grid + seq_len(n_grid)
        state_jacobian[, state_columns[[component]] - n_parameters] <-
          transform_state_jacobian(weak_value$Ju[, columns, drop = FALSE], component)
      }
      parameter_jacobian <- sweep(weak_value$Jp, 2, coordinates$pscale, "*")
      ode_jacobian <- sqrt(lambda) * weights$W %*%
        cbind(parameter_jacobian, state_jacobian)
      residual_jacobian <- rbind(linear_jacobian, ode_jacobian)
      gradient <- as.vector(crossprod(residual_jacobian, residual))
    }
    list(
      r = residual, J = residual_jacobian, value = sum(contributions), gradient = gradient,
      p = p, U = U, contributions = contributions,
      prior_contributions = setNames(colSums(prior_residual^2) / 2, prior_components),
      raw = weak_value$r
    )
  }
}

# Projected Levenberg-Marquardt for trajectory and ODE least squares.
# The augmented QR solve avoids squaring the Jacobian's condition number.
optimize_joint_gp <- function(par, evaluate, control, lower = rep(-Inf, length(par)),
                              upper = rep(Inf, length(par)), active = seq_along(par),
                              maxit = control$maxit) {
  par <- pmax(lower, pmin(upper, par))
  evaluations <- c(objective = 0L, jacobian = 0L)
  evaluate_active <- function(x, jacobian = TRUE) {
    full <- par
    full[active] <- x
    evaluations["objective"] <<- evaluations["objective"] + 1L
    if (jacobian) evaluations["jacobian"] <<- evaluations["jacobian"] + 1L
    value <- evaluate(full, jacobian)
    if (jacobian) value$J <- value$J[, active, drop = FALSE]
    value
  }
  x <- par[active]
  lower <- lower[active]
  upper <- upper[active]
  current <- evaluate_active(x)
  if (any(!is.finite(current$r)) || any(!is.finite(current$J))) {
    stop("Nonfinite initial LM residual or Jacobian.")
  }
  value <- sum(current$r^2) / 2
  damping <- 1e-5
  converged <- FALSE
  reason <- "iteration limit"
  for (iteration in seq_len(maxit)) {
    gradient <- as.vector(crossprod(current$J, current$r))
    scale <- sqrt(colSums(current$J^2))
    scale[scale == 0] <- 1 # An exactly inactive column has zero gradient and step.
    blocked <- (x <= lower & gradient > 0) | (x >= upper & gradient < 0)
    projected <- gradient
    projected[blocked] <- 0
    if (max(abs(projected) / scale) <= control$gtol) {
      converged <- TRUE
      reason <- "scaled gradient"
      break
    }
    free <- which(!blocked)
    scaled_jacobian <- sweep(current$J[, free, drop = FALSE], 2, scale[free], "/")
    augmented <- rbind(scaled_jacobian, diag(sqrt(damping), length(free)))
    target <- c(-current$r, rep(0, length(free)))
    step <- numeric(length(x))
    step[free] <- as.vector(qr.solve(augmented, target, tol = 1e-12)) / scale[free]
    candidate <- pmax(lower, pmin(upper, x + step))
    step <- candidate - x
    predicted <- -sum(gradient * step) - sum((current$J %*% step)^2) / 2
    trial <- evaluate_active(candidate, FALSE)
    trial_value <- if (all(is.finite(trial$r))) sum(trial$r^2) / 2 else Inf
    gain <- if (predicted > 0) (value - trial_value) / predicted else -Inf
    if (is.finite(gain) && gain > 1e-4 && trial_value < value) {
      small_value <- value - trial_value <= control$ftol * max(1, value)
      small_step <- sqrt(sum((scale * step)^2)) <=
        control$xtol * (control$xtol + sqrt(sum((scale * x)^2)))
      x <- candidate
      value <- trial_value
      current <- evaluate_active(x)
      if (any(!is.finite(current$J))) stop("Nonfinite LM Jacobian after an accepted step.")
      damping <- max(1e-12, damping * max(1 / 3, 1 - (2 * gain - 1)^3))
      if (small_value || small_step) {
        converged <- TRUE
        reason <- if (small_step) "step tolerance" else "objective tolerance"
        break
      }
    } else {
      damping <- min(1e16, damping * 5)
    }
    if (damping >= 1e16) {
      reason <- "damping limit"
      break
    }
  }
  par[active] <- x
  list(
    par = par, objective = value, method = "lm", converged = converged,
    code = if (converged) 0L else 1L, reason = reason,
    iterations = iteration, evaluations = evaluations
  )
}

# Initialize parameters with the state fixed at its GP mean. The first pass
# uses unweighted interior residuals; later passes use the test-Gram metric.
initialize_ode_parameters <- function(U, p0, weak, model, control, lower, upper) {
  rows <- weak_component_rows(seq_len(weak$ni), weak$K, model$D)
  unweighted <- diag(length(rows))
  weighted <- build_covariance_whitener(build_ode_metric(weak, model$D, rows), control$covariance_tol)
  p <- p0
  for (iteration in 0:control$init_gls) {
    whitener <- if (iteration == 0) unweighted else weighted
    evaluate <- function(p, jacobian = TRUE) {
      weak_value <- evaluate_weak_residual(U, p, weak, model, jacobian)
      residual <- as.vector(whitener %*% weak_value$r[rows])
      jacobian <- if (jacobian) whitener %*% weak_value$Jp[rows, , drop = FALSE]
      list(
        r = residual, J = jacobian,
        value = sum(residual^2) / 2,
        gradient = if (!is.null(jacobian)) as.vector(crossprod(jacobian, residual))
      )
    }
    next_parameters <- optimize_joint_gp(p, evaluate, control, lower, upper)$par
    change <- max(abs(next_parameters - p) / pmax(1, abs(p)))
    p <- next_parameters
    if (iteration > 0 && change < control$xtol) break
  }
  p
}

# ---- Observed formulation --------------------------------------------------
# The ODE precision is multiplied by n * mean(rms^2 / noise^2) so it keeps
# pace with the likelihood as data grow. There are no GP prior rows: with the
# ODE enforced they act as a parameter prior that is biased when the GP
# misfits the trajectory.
fit_observed_joint_gp <- function(f, Y, tt, p0, noise_sd, control, lower, upper) {
  n_components <- ncol(Y)
  if (control$verbose) message("Fitting ", n_components, " component GP(s).")
  fits <- lapply(seq_len(n_components), function(component) {
    fit_component_gp(tt, Y[, component], noise_sd[component], control)
  })
  noise <- sqrt(vapply(fits, `[[`, numeric(1), "noise2"))
  grid <- build_joint_gp_grid(tt, control)
  state <- prepare_gp_state(fits, grid, control)
  H <- build_observation_matrix(grid, tt)
  model <- build_symbolic_ode_model(f, n_components, length(p0))
  design <- build_weak_test_pool(grid, control)
  weak <- build_gauss_weak_operator(grid, design, fits, state$L, control, build_test_basis(grid, design, control))

  pilot <- initialize_ode_parameters(state$mean, p0, weak, model, control, lower, upper)
  metric <- build_ode_metric(weak, n_components)
  signal_to_noise <- mean(colMeans(Y^2) / noise^2)
  precision <- calculate_ode_scaling(Y, tt, control)$precision * length(tt) * signal_to_noise
  weights <- build_ode_weights(metric, weak, n_components, control, precision)

  coordinates <- state$coordinates
  coordinates$pscale <- rep(1, model$J)
  objective <- make_joint_gp_objective(Y, H, weak, model, weights, noise, control$lambda, coordinates)
  bounds <- build_optimizer_bounds(coordinates, lower, upper, state$bounds)
  start <- pack_joint_gp_coordinates(pilot, state$mean, coordinates)
  result <- optimize_joint_gp(start, objective, control, bounds$lower, bounds$upper)
  point <- unpack_joint_gp_coordinates(result$par, coordinates)
  final <- objective(result$par, FALSE)
  if (!result$converged) warning(result$reason, " converged=FALSE.", call. = FALSE)
  structure(list(
    phat = point$p, U_hat = point$U, tt = grid, tt_obs = tt, Y = Y,
    observed = seq_len(n_components), latent = integer(), gp = fits, noise_sd = noise,
    objective = result$objective, contributions = final$contributions,
    converged = result$converged, convergence_reason = result$reason,
    optimizer = result[setdiff(names(result), "par")],
    problem = list(
      evaluate = objective, theta = result$par, coordinates = coordinates, H = H,
      weak = weak, model = model, weights = weights, control = control
    )
  ), class = "jointgp")
}

# ---- Latent formulation ----------------------------------------------------
# One entirely unobserved component, interpolated in the weak form by a fixed
# GP. Observed components keep their data and GP prior.

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

# Observed columns are interpolated to the grid; the latent column starts at
# their row mean (or the largest-RMS column if the mean cancels).
initialize_latent_state <- function(U, tt, grid, observed, latent, component_scale = NULL) {
  initial_state <- matrix(0, length(grid), ncol(U))
  for (component in observed) {
    initial_state[, component] <- stats::approx(tt, U[, component], xout = grid, rule = 2)$y
  }
  observed_state <- initial_state[, observed, drop = FALSE]
  if (!is.null(component_scale)) observed_state <- sweep(observed_state, 2, component_scale[observed], "/")
  curve <- rowMeans(observed_state)
  observed_rms <- sqrt(colMeans(observed_state^2))
  if (sqrt(mean(curve^2)) <= sqrt(.Machine$double.eps) * max(observed_rms)) {
    curve <- observed_state[, which.max(observed_rms)]
  }
  if (!is.null(component_scale)) curve <- curve * component_scale[latent]
  initial_state[, latent] <- curve
  initial_state
}

predict_gp_mean <- function(fit, tt) {
  as.vector(fit$mean + gp_cross_covariance(fit, tt, fit$tt) %*% fit$alpha)
}

# WSINDy time scale (Messenger & Bortz, arXiv:2007.02848, eq. 4.7): the
# largest test has unit derivative sup-norm in t̃ = γ_t t, so γ_t = ‖φ′‖∞ / r_max.
calculate_test_time_scale <- function(design, eta) {
  max(abs(evaluate_bump(seq(0, 1, length.out = 1001L), 1L, eta))) / max(design$radii)
}

# Observed components use their fitted GPs. The latent GP is a fixed
# stationary Matern with unit radius in gamma_t time (equal radius bounds),
# amplitude units[latent] and the pilot mean.
build_latent_fits <- function(problem) {
  fits <- vector("list", problem$n_components)
  fits[problem$observed] <- problem$observed_fits
  fits[[problem$latent]] <- list(
    origin = min(problem$tt), span = 1 / problem$time_scale,
    tau2 = problem$units[problem$latent]^2, radius_coef = 0,
    radius_bounds = c(1, 1), mean = problem$latent_mean
  )
  fits
}

# Observed blocks are optimized in fixed GP posterior coordinates given the
# data. The objective still uses the original likelihood and prior.
build_posterior_coordinates <- function(problem) {
  n_grid <- length(problem$grid)
  prior <- problem$observed_prior
  coordinates <- prior
  for (k in seq_along(problem$observed)) {
    component <- problem$observed[k]
    L <- prior$L[[component]]
    A <- problem$H %*% L / problem$noise[k]
    R <- chol(diag(n_grid) + crossprod(A))
    R_inv <- backsolve(R, diag(n_grid))
    coordinates$L[[component]] <- t(chol(tcrossprod(L %*% R_inv)))
    coordinates$mu[, component] <- prior$mu[, component] + L %*% solve_cholesky_system(
      R,
      crossprod(A, (problem$Y[, k] - problem$H %*% prior$mu[, component]) / problem$noise[k])
    )
  }
  coordinates
}

# Fixed parameter scales for the latent optimizer coordinates.
resolve_parameter_scale <- function(p0, control) {
  scale <- if (is.null(control$parameter_scale)) ifelse(p0 == 0, 1, abs(p0)) else control$parameter_scale
  if (length(scale) != length(p0)) stop("parameter_scale must have one value per parameter.")
  scale
}

# Data, observed GPs, grid, pilot state, tests, component units, and the fixed
# weak operator with its weights.
prepare_latent_problem <- function(f, U, tt, p0, lower, upper, noise_sd, control, latent_initial) {
  components <- identify_latent_component(U)
  observed <- components$observed
  latent <- components$latent
  n_components <- ncol(U)
  Y <- U[, observed, drop = FALSE]
  observed_fits <- lapply(seq_along(observed), function(k) {
    fit_component_gp(tt, Y[, k], noise_sd[k], control)
  })
  grid <- build_joint_gp_grid(tt, control)
  component_scale <- control$ode_component_scale
  if (!is.null(component_scale)) {
    if (!length(component_scale) %in% c(1L, n_components)) {
      stop("ode_component_scale must be scalar or have one value per component.")
    }
    component_scale <- rep_len(component_scale, n_components)
  }
  scale_pilot <- initialize_latent_state(U, tt, grid, observed, latent, component_scale)
  initial_state <- scale_pilot
  if (!is.null(latent_initial)) {
    initial_state[, latent] <- stats::approx(tt, latent_initial, xout = grid, rule = 2)$y
  }
  bounds <- attr(control, "state_bounds")
  # Bounded observed starts use the GP mean: interpolated raw noise can leave
  # the model's domain.
  for (k in seq_along(observed)) {
    if (bounds$bounded[observed[k]]) {
      initial_state[, observed[k]] <- predict_gp_mean(observed_fits[[k]], grid)
    }
  }
  initial_state <- clip_state_to_bounds(initial_state, bounds)
  design <- build_weak_test_pool(grid, control)
  maps <- build_test_basis(grid, design, control)
  # Units are fixed from the observed RMS and the default latent initializer,
  # independently of latent_initial.
  units <- sqrt(colMeans(scale_pilot^2))
  units[observed] <- sqrt(colMeans(Y^2))
  if (any(!is.finite(units) | units <= 0)) stop("Unidentified pilot component scale.")
  scaling <- calculate_ode_scaling(matrix(rep(units, each = length(tt)), length(tt)), tt, control)
  problem <- list(
    Y = Y, tt = tt, grid = grid, H = build_observation_matrix(grid, tt), observed_fits = observed_fits,
    noise = sqrt(vapply(observed_fits, `[[`, numeric(1), "noise2")),
    model = build_symbolic_ode_model(f, n_components, length(p0)),
    units = scaling$rms, pilot = initial_state, latent_mean = mean(scale_pilot[, latent]),
    time_scale = calculate_test_time_scale(design, control$bump_eta),
    n_components = n_components, observed = observed, latent = latent,
    control = control, p0 = p0, lower = lower, upper = upper
  )

  problem$fits <- build_latent_fits(problem)
  problem$factors <- lapply(problem$fits, function(fit) {
    t(factor_gp_covariance(gp_cross_covariance(fit, grid, grid), control$gp_jitter)$R)
  })
  weak <- build_gauss_weak_operator(grid, design, problem$fits, problem$factors, control, maps, latent)
  problem$parameter_scale <- resolve_parameter_scale(p0, control)
  metric <- build_ode_metric(weak, n_components)
  problem$weak <- weak
  problem$weights <- build_ode_weights(metric, weak, n_components, control, scaling$precision)

  observed_prior <- list(mu = matrix(0, length(grid), n_components), L = vector("list", n_components))
  for (k in seq_along(observed)) {
    observed_prior$mu[, observed[k]] <- observed_fits[[k]]$mean
    observed_prior$L[[observed[k]]] <- problem$factors[[observed[k]]]
  }
  problem$observed_prior <- observed_prior
  problem$posterior_coordinates <- build_posterior_coordinates(problem)
  problem
}

# The latent block is optimized in scaled physical values; only observed
# components have GP prior rows.
make_latent_gp_objective <- function(problem) {
  latent <- problem$latent
  coordinates <- problem$posterior_coordinates
  coordinates$mu[, latent] <- 0
  coordinates$L[[latent]] <- diag(problem$units[latent], length(problem$grid))
  coordinates$order <- c(problem$observed, latent)
  coordinates$pscale <- problem$parameter_scale
  state_bounds <- attr(problem$control, "state_bounds")
  coordinates <- apply_state_bounds_to_coordinates(coordinates, state_bounds, problem$units)
  evaluate <- make_joint_gp_objective(
    problem$Y, problem$H, problem$weak, problem$model, problem$weights,
    problem$noise, problem$control$lambda, coordinates, problem$observed_prior, problem$observed, problem$observed
  )
  objective <- function(par, jacobian = TRUE) {
    value <- evaluate(par, jacobian)
    contributions <- value$contributions
    value$contributions <- c(
      data = unname(contributions["data"]),
      observed_gp = sum(value$prior_contributions[as.character(problem$observed)]),
      ode = sum(contributions[c("interior", "boundary")])
    )
    value
  }
  attr(objective, "coordinates") <- coordinates
  attr(objective, "bounds") <- build_optimizer_bounds(coordinates, problem$lower, problem$upper, state_bounds)
  objective
}

fit_latent_joint_gp <- function(f, U, tt, p0, noise_sd, control, latent_initial, lower, upper) {
  problem <- prepare_latent_problem(f, U, tt, p0, lower, upper, noise_sd, control, latent_initial)
  n_parameters <- length(p0)
  n_grid <- length(problem$grid)
  objective <- make_latent_gp_objective(problem)
  coordinates <- attr(objective, "coordinates")
  bounds <- attr(objective, "bounds")
  par <- pack_joint_gp_coordinates(p0, problem$pilot, coordinates)
  par <- pmax(bounds$lower, pmin(bounds$upper, par))
  # Warm start on the parameters and the latent block, with the observed blocks
  # held at their start.
  latent_block <- n_parameters + n_grid * length(problem$observed) + seq_len(n_grid)
  par <- optimize_joint_gp(par, objective, control,
    maxit = min(75L, control$maxit),
    active = c(seq_len(n_parameters), latent_block),
    lower = bounds$lower, upper = bounds$upper
  )$par
  result <- optimize_joint_gp(par, objective, control, bounds$lower, bounds$upper)
  final <- objective(result$par, FALSE)
  if (!is.finite(final$value) || any(!is.finite(c(final$p, final$U)))) stop("No finite latent estimate was produced.")
  if (!result$converged) warning(result$reason, " converged=FALSE.", call. = FALSE)
  structure(list(
    phat = final$p, U_hat = final$U, tt = problem$grid, tt_obs = tt, Y = problem$Y,
    observed = problem$observed, latent = problem$latent, gp = problem$observed_fits, noise_sd = problem$noise,
    objective = final$value, contributions = final$contributions,
    converged = result$converged, convergence_reason = result$reason,
    optimizer = result[setdiff(names(result), "par")],
    problem = list(
      evaluate = objective, theta = result$par, coordinates = coordinates, H = problem$H,
      weak = problem$weak, model = problem$model, weights = problem$weights, control = control
    )
  ), class = c("jointgp", "wendygp_latent"))
}
