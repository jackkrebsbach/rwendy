# Multiscale weak test selection by SVD or parameter sensitivity.

weak_component_rows <- function(test_indices, n_tests, n_components) {
  unlist(lapply(seq_len(n_components), function(component) (component - 1L) * n_tests + test_indices), use.names = FALSE)
}

build_weak_test_pool <- function(grid, control, radii = NULL) {
  span <- diff(range(grid))
  start <- min(grid)
  end <- max(grid)
  if (is.null(radii)) radii <- control$weak_design_radii
  if (is.null(radii)) {
    # Halve the radius down to one grid spacing, the smallest resolved support.
    spacing <- span / max(1L, length(grid) - 1L)
    radii <- .4 * span / 2^(0:max(0L, floor(log2(.4 * span / spacing))))
  }
  radii <- sort(unique(radii))
  if (any(radii >= span / 2)) stop("weak_design_radii must be smaller than half the span.")
  interior <- do.call(rbind, lapply(radii, function(radius) {
    if (is.null(control$weak_design_centers)) {
      # No subsampling or added off-grid endpoints. The small tolerance only
      # admits exact support endpoints affected by floating-point roundoff.
      candidates <- grid
      tol <- 32 * .Machine$double.eps * max(1, abs(start), abs(end), span)
      centers <- candidates[candidates >= start + radius - tol & candidates <= end - radius + tol]
    } else {
      n_centers <- max(3L, min(control$weak_design_centers, ceiling((span - 2 * radius) / (radius / 4)) + 1L))
      centers <- seq(start + radius, end - radius, length.out = n_centers)
    }
    data.frame(center = centers, radius = rep(radius, length(centers)))
  }))
  if (!nrow(interior)) stop("No admissible interior grid centers; refine the working grid.")
  boundary_radii <- if (is.null(control$bl_radii)) max(radii) else control$bl_radii
  if (any(boundary_radii >= span / 2)) stop("bl_radii must be smaller than half the span.")
  boundary <- if (control$include_bl) {
    do.call(rbind, lapply(boundary_radii, function(radius) {
      offsets <- seq(0, .4 * radius, length.out = control$bl_count)
      data.frame(center = c(start + offsets, end - offsets), radius = radius)
    }))
  } else {
    data.frame(center = numeric(), radius = numeric())
  }
  list(
    interior = interior, boundary = boundary, radii = radii, boundary_radii = boundary_radii,
    radius_selection = list(
      method = "geometric",
      center_placement = if (is.null(control$weak_design_centers)) "grid" else "capped",
      grid_points = length(grid)
    )
  )
}

# Compare the complete weak residual AND parameter sensitivities on nested
# quadratures. Normalize by integral magnitudes, not measurement noise or a
# residual that should be zero at a solution. Screening is local at the pilot.
measure_weak_test_errors <- function(weak, fine, U, fine_U, p, model, pscale) {
  coarse_value <- evaluate_weak_residual(U, p, weak, model)
  fine_value <- evaluate_weak_residual(fine_U, p, fine, model)
  rhs <- model$jet[[1]](rbind(matrix(p, model$J, nrow(U)), t(U), weak$tt))
  magnitude <- abs(weak$V) %*% abs(rhs) + abs(weak$Vp) %*% abs(U)
  delta <- matrix(coarse_value$r - fine_value$r, weak$K, model$D)
  sensitivity_difference <- sweep(coarse_value$Jp - fine_value$Jp, 2, pscale, "*")
  sensitivity <- matrix(sqrt(rowSums(sensitivity_difference^2)), weak$K, model$D)
  for (component in seq_len(model$D)) {
    magnitude[, component] <- pmax(
      magnitude[, component],
      .Machine$double.eps * max(magnitude[, component]), .Machine$double.xmin
    )
  }
  list(
    base = coarse_value, delta = as.vector(delta),
    relative = apply(pmax(abs(delta), sensitivity) / magnitude, 1, max)
  )
}

# Propagate the state covariance conditional on the observations through Ju.
propagate_posterior_state_covariance <- function(Ju, state, H, noise) {
  n_grid <- nrow(state$mu)
  covariance <- matrix(0, nrow(Ju), nrow(Ju))
  for (component in seq_along(state$L)) {
    prior_factor <- state$L[[component]]
    observation_factor <- H %*% prior_factor / noise[component]
    precision <- chol(diag(n_grid) + crossprod(observation_factor))
    columns <- (component - 1L) * n_grid + seq_len(n_grid)
    residual_factor <- Ju[, columns, drop = FALSE] %*% prior_factor
    factor <- forwardsolve(t(precision), t(residual_factor))
    covariance <- covariance + crossprod(factor)
  }
  covariance
}

prepare_weak_test_selection <- function(tt, grid, fits, state, H, p, model, control, pscale,
                                        noise) {
  design <- build_weak_test_pool(grid, control)
  n_interior <- nrow(design$interior)
  n_boundary <- nrow(design$boundary)
  raw <- build_grid_weak_operator(
    grid, design, control,
    list(interior = diag(n_interior), boundary = diag(n_boundary))
  )
  fine_times <- seq(min(grid), max(grid), length.out = 2L * length(grid) - 1L)
  fine_mean <- do.call(cbind, lapply(fits, function(g) predict_component_gp(g, fine_times)$mean))
  fine <- build_grid_weak_operator(fine_times, design, control, raw$maps)
  check <- measure_weak_test_errors(raw, fine, state$mean, fine_mean, p, model, pscale)
  keep_raw <- which(check$relative[seq_len(n_interior)] <= control$weak_design_quad_tol)
  fallback <- !length(keep_raw)
  if (fallback) keep_raw <- which(design$interior$radius == max(design$radii))
  basis <- orthonormalize_weak_tests(evaluate_weak_test_rows(
    grid, design$interior[keep_raw, , drop = FALSE],
    0L, control$bump_eta
  ), control$basis_tol, spectrum = TRUE)
  basis_map <- basis$map
  map <- matrix(0, nrow(basis_map), n_interior)
  map[, keep_raw] <- basis_map
  boundary_map <- orthonormalize_weak_tests(
    evaluate_weak_test_rows(grid, design$boundary, 0L, control$bump_eta), control$basis_tol
  )
  weak <- build_grid_weak_operator(grid, design, control, list(interior = map, boundary = boundary_map))
  fine <- build_grid_weak_operator(fine_times, design, control, weak$maps)
  check_modes <- measure_weak_test_errors(weak, fine, state$mean, fine_mean, p, model, pscale)
  weak_value <- check_modes$base
  ode_covariance <- build_ode_metric(weak, model$D, control, weak_value$Ju, state$Sigma)
  sd <- sqrt(pmax(diag(ode_covariance), .Machine$double.xmin))
  standard_error <- apply(matrix(abs(check_modes$delta) / sd, weak$K, model$D), 1, max)
  keep_modes <- which(check_modes$relative[seq_len(weak$ni)] <= control$weak_design_quad_tol &
    standard_error[seq_len(weak$ni)] <= control$weak_grid_tol)
  if (!length(keep_modes)) {
    keep_modes <- 1L
    fallback <- TRUE
  }
  physical <- c(keep_modes, weak$ni + seq_len(weak$nb))
  rows <- weak_component_rows(physical, weak$K, model$D)
  maps <- list(interior = weak$maps$interior[keep_modes, , drop = FALSE], boundary = weak$maps$boundary)
  selected <- build_grid_weak_operator(grid, design, control, maps)
  Ju <- weak_value$Ju[rows, , drop = FALSE]
  singular_values <- basis$singular_values[basis$retained[keep_modes]]
  singular_total <- sum(basis$singular_values)
  information_available <- sum(singular_values) / singular_total
  information_passed <- information_available + 32 * .Machine$double.eps >= control$weak_design_info
  list(
    weak = selected, G = weak_value$Jp[rows, , drop = FALSE], Ju = Ju,
    Omega = ode_covariance[rows, rows, drop = FALSE],
    ode_precision = calculate_ode_scaling(do.call(cbind, lapply(fits, `[[`, "y")), tt, control)$precision,
    singular_values = singular_values, singular_total = singular_total,
    include_gp_prior = !isFALSE(state$include_gp_prior),
    nuisance = if (control$weak_radius_method == "sensitivity" && state$include_gp_prior) {
      propagate_posterior_state_covariance(Ju, state, H, noise)
    } else {
      NULL
    },
    data_jacobian = if (control$weak_radius_method == "sensitivity" && !state$include_gp_prior) {
      kronecker(diag(1 / noise, length(noise)), H)
    } else {
      NULL
    },
    screening = list(
      passed = !fallback && (!is.null(control$weak_design_budget) || information_passed),
      center_placement = design$radius_selection$center_placement,
      grid_points = length(grid), raw_count = n_interior, raw_retained = length(keep_raw),
      mode_count = weak$ni, mode_retained = length(keep_modes),
      raw_relative_error = check$relative[seq_len(n_interior)],
      mode_relative_error = check_modes$relative[seq_len(weak$ni)],
      mode_standardized_error = standard_error[seq_len(weak$ni)],
      retained_raw = keep_raw, retained_modes = keep_modes,
      singular_values = basis$singular_values, information_available = information_available,
      information_target = control$weak_design_info, information_passed = information_passed,
      tolerance = control$weak_design_quad_tol, fallback = fallback
    )
  )
}

# MSG information-number convention: cumulative singular values, not their
# squares. Never renormalize away modes lost to numerical/accuracy screening.
select_singular_value_count <- function(values, total, target) {
  if (!length(values) || !is.finite(total) || total <= 0) {
    stop("No finite positive singular-value information is available.")
  }
  fraction <- cumsum(values) / total
  reached <- which(fraction + 32 * .Machine$double.eps >= target)
  list(
    count = if (length(reached)) reached[1L] else length(values),
    fractions = fraction, target_met = length(reached) > 0L
  )
}

# Profile completely free grid-state increments out of the stacked data/ODE
# Jacobian. This is a generalized Schur complement even when the data-only
# state Hessian is singular. It uses no GP covariance or ranking ridge.
profile_state_jacobian <- function(parameter_jacobian, state_jacobian) {
  scales <- pmax(sqrt(colSums(state_jacobian^2)), .Machine$double.xmin)
  decomposition <- svd(sweep(state_jacobian, 2, scales, "/"), nu = min(dim(state_jacobian)), nv = 0)
  tol <- max(dim(state_jacobian)) * .Machine$double.eps
  keep <- which(decomposition$d > max(decomposition$d) * tol)
  if (length(keep) == nrow(state_jacobian)) {
    residual <- parameter_jacobian * 0
  } else {
    state_basis <- decomposition$u[, keep, drop = FALSE]
    residual <- parameter_jacobian - state_basis %*% crossprod(state_basis, parameter_jacobian)
  }
  list(
    matrix = crossprod(residual), state_rank = length(keep),
    state_nullity = ncol(state_jacobian) - length(keep), state_rank_tolerance = tol
  )
}

# EXACT Schur complement of the local joint Gauss-Newton matrix (up to the
# existing covariance rank projection), computed in residual rather than state
# dimension. Parameter scales are fixed from user inputs, never truth.
measure_parameter_information <- function(space, selected, model, lambda, pscale, control) {
  weak <- space$weak
  physical <- c(sort(selected), weak$ni + seq_len(weak$nb))
  rows <- weak_component_rows(physical, weak$K, model$D)
  layout <- list(K = length(physical), ni = length(selected), nb = weak$nb)
  weights <- build_ode_weights(space$Omega[rows, rows, drop = FALSE], layout, model$D, control)
  if (!is.null(space$ode_precision)) weights <- scale_weak_weights(weights, space$ode_precision, layout$K)
  parameter_jacobian <- sqrt(lambda) * sweep(weights$W %*% space$G[rows, , drop = FALSE], 2, pscale, "*")
  free <- NULL
  if (isFALSE(space$include_gp_prior)) {
    state_jacobian <- rbind(space$data_jacobian, sqrt(lambda) * weights$W %*% space$Ju[rows, , drop = FALSE])
    free <- profile_state_jacobian(
      rbind(matrix(0, nrow(space$data_jacobian), ncol(parameter_jacobian)), parameter_jacobian),
      state_jacobian
    )
    information <- free$matrix
  } else {
    marginal_covariance <- diag(nrow(weights$W)) + lambda * weights$W %*%
      space$nuisance[rows, rows, drop = FALSE] %*% t(weights$W)
    covariance_factor <- chol((marginal_covariance + t(marginal_covariance)) / 2)
    information <- crossprod(forwardsolve(t(covariance_factor), parameter_jacobian))
  }
  c(
    list(
      matrix = information, eigenvalues = pmax(0, eigen(information, symmetric = TRUE, only.values = TRUE)$values),
      covariance_rank = nrow(weights$W), include_gp_prior = !isFALSE(space$include_gp_prior)
    ),
    if (!is.null(free)) free[setdiff(names(free), "matrix")]
  )
}

select_weak_tests <- function(space, model, control, pscale) {
  count <- space$weak$ni
  spectral <- select_singular_value_count(
    space$singular_values, space$singular_total,
    control$weak_design_info
  )
  budget <- if (is.null(control$weak_design_budget)) {
    spectral$count
  } else {
    min(count, control$weak_design_budget)
  }
  selected <- seq_len(budget)
  information <- NULL
  if (control$weak_radius_method == "sensitivity") {
    full <- measure_parameter_information(
      space, seq_len(count), model,
      control$lambda, pscale, control
    )
    # A common floor lets us compare singular candidate information matrices.
    floor <- max(max(full$eigenvalues) * 1e-8, .Machine$double.eps)
    score <- function(candidate) {
      info <- measure_parameter_information(
        space, candidate, model,
        control$lambda, pscale, control
      )
      -sum(1 / (info$eigenvalues + floor))
    }
    selected <- seq_len(min(budget, control$weak_design_coverage))
    while (length(selected) < budget) {
      remaining <- setdiff(seq_len(count), selected)
      scores <- vapply(remaining, function(j) score(c(selected, j)), numeric(1))
      selected <- c(selected, remaining[which.max(scores)])
    }
    selected <- sort(selected)
    information <- measure_parameter_information(
      space, selected, model,
      control$lambda, pscale, control
    )
  }
  fraction <- sum(space$singular_values[selected]) / space$singular_total
  diagnostics <- list(
    method = control$weak_radius_method, selected = selected,
    budget = budget, information_fraction = fraction,
    information_target = control$weak_design_info,
    information_target_met = fraction + 32 * .Machine$double.eps >= control$weak_design_info,
    singular_values = space$singular_values, singular_total = space$singular_total,
    information = information, screening = space$screening
  )
  maps <- list(
    interior = space$weak$maps$interior[selected, , drop = FALSE],
    boundary = space$weak$maps$boundary
  )
  design <- space$weak$design
  design$radius_selection <- diagnostics
  weak <- build_grid_weak_operator(space$weak$tt, design, control, maps)
  physical <- c(selected, space$weak$ni + seq_len(space$weak$nb))
  rows <- weak_component_rows(physical, space$weak$K, model$D)
  list(weak = weak, Omega = space$Omega[rows, rows, drop = FALSE], diagnostics = diagnostics)
}

# Prepare and fit the requested weak space on one fixed working grid.
solve_with_selected_weak_tests <- function(Y, tt, p0, fits, model, control, lower, upper) {
  pscale <- control$weak_design_scale
  if (is.null(pscale)) pscale <- pmax(abs(p0), 1)
  if (length(pscale) != model$J) stop("weak_design_scale must have one entry per parameter.")
  grid <- build_joint_gp_grid(tt, control)
  noise <- sqrt(vapply(fits, `[[`, numeric(1), "noise2"))
  state <- prepare_gp_state(fits, grid, control)
  H <- build_observation_matrix(grid, tt)
  pilot_design <- build_weak_test_pool(grid, control, diff(range(tt)) * c(.25, .4))
  pilot_weak <- build_grid_weak_operator(grid, pilot_design, control)
  pilot <- initialize_ode_parameters(
    state$mean, p0, pilot_weak, model, state$Sigma,
    control, lower, upper
  )
  space <- prepare_weak_test_selection(
    tt, grid, fits, state, H, pilot, model,
    control, pscale, noise
  )
  spec <- select_weak_tests(space, model, control, pscale)
  weights <- prepare_ode_metric(spec$Omega, spec$weak, Y, tt, control)$weights
  initial_error <- check_initial_grid_refinement(state, pilot, spec$weak, model, fits, weights, control)
  interpolation <- vapply(fits, measure_gp_interpolation_error, numeric(1),
    tt = grid, obs = tt, H = H
  )
  fit <- fit_observed_joint_gp(
    Y, tt, fits, model, control, lower, upper,
    state, H, pilot, spec, initial_error
  )
  fit$diagnostics$grid_passed <- fit$diagnostics$grid_passed && max(interpolation) <= control$grid_tol
  fit$diagnostics$initial_interpolation <- interpolation
  fit
}
