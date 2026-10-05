# Integrate the state extension with Gauss quadrature on each grid interval.
# Quadrature nodes are not additional optimization variables.
build_gauss_rule <- function(grid, order) {
  indices <- seq_len(order - 1L)
  off_diagonal <- indices / sqrt(4 * indices^2 - 1)
  jacobi_matrix <- matrix(0, order, order)
  jacobi_matrix[cbind(indices, indices + 1L)] <- off_diagonal
  jacobi_matrix <- jacobi_matrix + t(jacobi_matrix)
  decomposition <- eigen(jacobi_matrix, symmetric = TRUE)
  sorted_indices <- order(decomposition$values)
  nodes <- decomposition$values[sorted_indices]
  weights <- 2 * decomposition$vectors[1, sorted_indices]^2 / sum(decomposition$vectors[1, sorted_indices]^2)
  interval_widths <- diff(grid)
  midpoints <- (head(grid, -1L) + tail(grid, -1L)) / 2
  list(
    tt = as.vector(t(outer(midpoints, rep(1, order)) + outer(interval_widths / 2, nodes))),
    weights = as.vector(t(outer(interval_widths / 2, weights))), order = order
  )
}

evaluate_gauss_test_rows <- function(tt, design, derivative, eta) {
  rbind(
    evaluate_weak_test_rows(tt, design$interior, derivative, eta),
    evaluate_weak_test_rows(tt, design$boundary, derivative, eta)
  )
}

# Tie the stencil width to the fitted GP length scale. A wider stencil can make
# high-order interpolation diverge. The working grid has uniform spacing.
choose_lagrange_stencil_size <- function(fit, grid, tt) {
  spacing <- mean(diff(grid))
  n_grid <- length(grid)
  # Each integration cell must carry one polynomial. Choosing the width at
  # each quadrature node introduces stencil switches inside a cell, hence
  # discontinuities that an unsplit Gauss rule resolves poorly.
  midpoints <- (head(grid, -1L) + tail(grid, -1L)) / 2
  radius <- evaluate_gp_radius((midpoints - fit$origin) / fit$span, fit$radius_coef, fit$radius_bounds)$r * fit$span
  width <- pmin(2L * pmax(2L, pmin(8L, as.integer(round(radius / spacing / 2)))), n_grid)
  width[findInterval(tt, grid, rightmost.closed = TRUE, all.inside = TRUE)]
}

# Centered-stencil Lagrange rows. Weights sum to one, so the constant mean
# offset the GP extension needs is a no-op here and the code shape is shared.
build_lagrange_interpolation <- function(grid, tt, stencil_sizes) {
  n_grid <- length(grid)
  interpolation <- matrix(0, length(tt), n_grid)
  interval <- findInterval(tt, grid, rightmost.closed = TRUE)
  for (i in seq_along(tt)) {
    width <- stencil_sizes[i]
    first <- min(max(interval[i] - width %/% 2L + 1L, 1L), n_grid - width + 1L)
    columns <- first + seq_len(width) - 1L
    nodes <- grid[columns]
    for (node in seq_len(width)) {
      others <- seq_len(width)[-node]
      interpolation[i, columns[node]] <- prod((tt[i] - nodes[others]) / (nodes[node] - nodes[others]))
    }
  }
  interpolation
}

# Per-component choice. Both are fixed representations in the weak operator;
# neither an invertible coordinate change nor an interpolant guarantees that
# an unobserved grid state is identifiable.
build_state_extensions <- function(fits, state, grid, tt, control) {
  if (!length(control$weak_extension) %in% c(1L, length(fits))) {
    stop("weak_extension must be scalar or have one entry per component.")
  }
  kind <- rep_len(control$weak_extension, length(fits))
  lapply(seq_along(fits), function(component) {
    if (identical(kind[component], "lagrange")) {
      build_lagrange_interpolation(grid, tt, choose_lagrange_stencil_size(fits[[component]], grid, tt))
    } else {
      t(solve_cholesky_system(t(state$L[[component]]), t(gp_cross_covariance(fits[[component]], tt, grid))))
    }
  })
}

build_gauss_weak_operator <- function(grid, design, fits, state, control, maps = NULL,
                                      order = control$weak_quad_order) {
  quad <- build_gauss_rule(grid, order)
  raw <- evaluate_gauss_test_rows(quad$tt, design, 0L, control$bump_eta)
  n_raw_interior <- nrow(design$interior)
  n_raw_boundary <- nrow(design$boundary)
  spectrum <- NULL
  if (is.null(maps)) {
    X <- sweep(raw, 2, sqrt(quad$weights), "*")
    spectrum <- orthonormalize_weak_tests(X[seq_len(n_raw_interior), , drop = FALSE], control$basis_tol, TRUE)
    size <- select_singular_value_count(
      spectrum$singular_values[spectrum$retained],
      sum(spectrum$singular_values), control$weak_design_info
    )
    count <- if (!is.null(control$weak_radii)) {
      nrow(spectrum$map)
    } else if (is.null(control$weak_design_budget)) {
      size$count
    } else {
      min(control$weak_design_budget, nrow(spectrum$map))
    }
    if (count < 1L) stop("Continuous weak test pool has no retained interior mode.")
    boundary_map <- if (n_raw_boundary) {
      orthonormalize_weak_tests(
        X[n_raw_interior + seq_len(n_raw_boundary), , drop = FALSE], control$basis_tol
      )
    } else {
      matrix(0, 0, 0)
    }
    maps <- list(interior = spectrum$map[seq_len(count), , drop = FALSE], boundary = boundary_map)
  }
  n_interior <- nrow(maps$interior)
  n_boundary <- nrow(maps$boundary)
  n_tests <- n_interior + n_boundary
  transform <- function(x) {
    rbind(
      maps$interior %*% x[seq_len(n_raw_interior), , drop = FALSE],
      maps$boundary %*% x[n_raw_interior + seq_len(n_raw_boundary), , drop = FALSE]
    )
  }
  phi <- transform(raw)
  prime <- transform(evaluate_gauss_test_rows(quad$tt, design, 1L, control$bump_eta))
  endpoints <- transform(evaluate_gauss_test_rows(range(grid), design, 0L, control$bump_eta))
  state_extensions <- build_state_extensions(fits, state, grid, quad$tt, control)
  boundary_extensions <- build_state_extensions(fits, state, grid, range(grid), control)
  weak <- list(
    integration = "gp_gauss", tt = grid, quad_tt = quad$tt, quad_weights = quad$weights,
    quad_order = order, ni = n_interior, nb = n_boundary, K = n_tests, design = design, maps = maps,
    V = sweep(phi, 2, quad$weights, "*"), Vp = sweep(prime, 2, quad$weights, "*"),
    B = cbind(-endpoints[, 1], endpoints[, 2]),
    gram = tcrossprod(sweep(phi, 2, sqrt(quad$weights), "*")),
    E = state_extensions, Eb = boundary_extensions, mean = vapply(fits, `[[`, numeric(1), "mean"), em_order = 0L,
    extension_kind = control$weak_extension
  )
  # These linear Jacobian terms stay fixed throughout the trajectory solve.
  weak$VpE <- lapply(state_extensions, function(E) weak$Vp %*% E)
  weak$BEb <- lapply(boundary_extensions, function(Eb) weak$B %*% Eb)
  if (!is.null(spectrum)) {
    singular_values <- spectrum$singular_values[spectrum$retained]
    total <- sum(spectrum$singular_values)
    fraction <- sum(head(singular_values, n_interior)) / total
    weak$selection <- list(
      method = if (is.null(control$weak_radii)) "svd" else "explicit", selected = seq_len(n_interior), budget = n_interior,
      size_rule = if (!is.null(control$weak_radii)) {
        "full_numerical_span"
      } else if (is.null(control$weak_design_budget)) {
        "singular_value_fraction"
      } else {
        "fixed_budget"
      },
      singular_values = singular_values, singular_total = total, information_fraction = fraction,
      information_target = control$weak_design_info,
      information_target_met = fraction + 32 * .Machine$double.eps >= control$weak_design_info,
      integration = "gp_gauss", screening = list(
        raw_count = n_raw_interior, mode_retained = n_interior,
        center_placement = design$radius_selection$center_placement,
        passed = fraction + 32 * .Machine$double.eps >= control$weak_design_info
      )
    )
  }
  weak
}

# Evaluate the weak residual and its Jacobian in physical state coordinates.
evaluate_gauss_weak_residual <- function(U, p, weak, model, jacobian = TRUE) {
  n_grid <- nrow(U)
  n_components <- model$D
  n_parameters <- model$J
  n_tests <- weak$K
  n_quadrature <- length(weak$quad_tt)
  quadrature_state <- matrix(0, n_quadrature, n_components)
  boundary_state <- matrix(0, 2L, n_components)
  for (component in seq_len(n_components)) {
    centered <- U[, component] - weak$mean[component]
    quadrature_state[, component] <- weak$mean[component] + weak$E[[component]] %*% centered
    boundary_state[, component] <- weak$mean[component] + weak$Eb[[component]] %*% centered
  }
  input <- rbind(matrix(p, n_parameters, n_quadrature), t(quadrature_state), weak$quad_tt)
  rhs <- model$jet[[1]](input)
  residual <- weak$V %*% rhs + weak$Vp %*% quadrature_state - weak$B %*% boundary_state
  parameter_jacobian <- NULL
  state_jacobian <- NULL
  if (jacobian) {
    rhs_jacobian <- array(model$jet_jac[[1]](input), c(n_quadrature, n_components, n_parameters + n_components))
    rhs_parameter_jacobian <- matrix(
      rhs_jacobian[, , seq_len(n_parameters), drop = FALSE],
      n_quadrature, n_components * n_parameters
    )
    parameter_jacobian <- matrix(
      weak$V %*% rhs_parameter_jacobian, n_tests * n_components, n_parameters
    )
    state_jacobian <- matrix(0, n_tests * n_components, n_grid * n_components)
    # Block (a, b) differentiates residual component a by state component b.
    for (a in seq_len(n_components)) {
      rows <- (a - 1L) * n_tests + seq_len(n_tests)
      for (b in seq_len(n_components)) {
        columns <- (b - 1L) * n_grid + seq_len(n_grid)
        # R recycles the derivative vector down each interpolation column.
        block <- weak$V %*% (weak$E[[b]] * rhs_jacobian[, a, n_parameters + b])
        if (a == b) {
          block <- block + weak$VpE[[b]] - weak$BEb[[b]]
        }
        state_jacobian[rows, columns] <- block
      }
    }
  }
  list(r = as.vector(residual), Jp = parameter_jacobian, Ju = state_jacobian)
}

solve_with_gauss_quadrature <- function(f, Y, tt, p0, fits, control, lower, upper) {
  grid <- build_joint_gp_grid(tt, control)
  state <- prepare_gp_state(fits, grid, control)
  H <- build_observation_matrix(grid, tt)
  n_components <- ncol(Y)
  n_parameters <- length(p0)
  model <- build_symbolic_ode_model(f, n_components, n_parameters, 0L)
  design <- if (is.null(control$weak_radii)) {
    build_weak_test_pool(grid, control)
  } else {
    build_weak_test_design(tt, fits, control)
  }
  weak <- build_gauss_weak_operator(grid, design, fits, state, control)
  weak$state_scale <- sqrt(vapply(fits, `[[`, numeric(1), "tau2"))
  weak$parameter_scale <- if (is.null(control$weak_design_scale)) pmax(abs(p0), 1) else control$weak_design_scale
  if (length(weak$parameter_scale) != n_parameters) stop("weak_design_scale must have one entry per parameter.")
  pilot <- initialize_ode_parameters(state$mean, p0, weak, model, state$Sigma, control, lower, upper)
  weak_value <- evaluate_weak_residual(state$mean, pilot, weak, model)
  ode_covariance <- build_ode_metric(weak, n_components, control, weak_value$Ju, state$Sigma)
  weights <- prepare_ode_metric(ode_covariance, weak, Y, tt, control)$weights
  initial <- check_gauss_quadrature_accuracy(state$mean, pilot, state, weak, model, fits, weights, H, control)
  weak_problem <- list(weak = weak, Omega = ode_covariance, diagnostics = weak$selection)
  weak_problem$diagnostics$screening$quadrature_passed <- initial$quadrature_passed
  weak_problem$diagnostics$screening$passed <- weak_problem$diagnostics$screening$passed && initial$quadrature_passed
  fit <- fit_observed_joint_gp(
    Y, tt, fits, model, control, lower, upper, state, H, pilot, weak_problem,
    initial$weak
  )
  interpolation_error <- vapply(fits, measure_gp_interpolation_error, numeric(1), tt = grid, obs = tt, H = H)
  fit$diagnostics$weak_grid_passed <- fit$diagnostics$final_grid$quadrature_passed
  fit$diagnostics$quadrature_passed <- fit$diagnostics$final_grid$quadrature_passed
  fit$diagnostics$grid_passed <- max(interpolation_error) <= control$grid_tol && fit$diagnostics$grid_passed
  fit$diagnostics$initial_quadrature <- initial
  fit$diagnostics$initial_interpolation <- interpolation_error
  fit$diagnostics$quadrature_points <- length(weak$quad_tt)
  fit$diagnostics$state_variables <- length(state$mu)
  if (!fit$diagnostics$grid_passed) {
    warning_text <- "Observation-operator error exceeds grid_tol on the fixed working grid."
    if (control$grid_action == "error") stop(warning_text) else warning(warning_text, call. = FALSE)
  }
  if (!fit$diagnostics$quadrature_passed) {
    warning(format_quadrature_warning(fit$diagnostics$final_grid, control), call. = FALSE)
  }
  if (!isTRUE(weak$selection$information_target_met)) {
    warning("Requested weak singular-value fraction was not reached; inspect radius_selection.", call. = FALSE)
  }
  fit
}
