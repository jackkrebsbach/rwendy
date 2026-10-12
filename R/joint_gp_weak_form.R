# Weak form on the working grid: grid and observation map, test functions,
# ODE model, Gauss-Legendre residual with state extensions, penalty weights.

# ---- Grid ------------------------------------------------------------------
# Piecewise-linear map from grid values to observation times.
build_observation_matrix <- function(grid, observation_times) {
  n_grid <- length(grid)
  observation_matrix <- matrix(0, length(observation_times), n_grid)
  left <- pmax(1L, pmin(n_grid - 1L, findInterval(observation_times, grid)))
  fraction <- (observation_times - grid[left]) / (grid[left + 1L] - grid[left])
  observation_matrix[cbind(seq_along(observation_times), left)] <- 1 - fraction
  observation_matrix[cbind(seq_along(observation_times), left + 1L)] <- fraction
  observation_matrix
}

# Smallest uniform grid with at least grid_min points and at least two grid
# intervals per observation interval. Regular observation times stay on it.
build_joint_gp_grid <- function(tt, control) {
  n_observations <- length(tt)
  spacing <- mean(diff(tt))
  regular <- max(abs(diff(tt) - spacing)) <= 1e-8 * spacing
  n_grid <- if (regular) {
    (n_observations - 1L) * max(2, ceiling((control$grid_min - 1L) / (n_observations - 1L))) + 1L
  } else {
    max(2L * n_observations - 1L, control$grid_min)
  }
  if (n_grid > control$grid_max) stop("grid_max cannot contain the fixed working grid.")
  seq(min(tt), max(tt), length.out = n_grid)
}

# ---- Test functions --------------------------------------------------------
# Unit-radius bump exp(eta - eta / (1 - x^2)) and its derivatives, scaled to
# peak 1 and exactly zero outside (-1, 1).
evaluate_bump <- function(x, order = 0L, eta = 9) {
  values <- x * 0
  inside <- abs(x) < 1 & is.finite(x)
  if (any(inside)) {
    bump <- function(t, r) phi(t, r, eta = eta)
    derivative <- test_function_derivative(bump, 1, 1, order)
    values[inside] <- exp(eta) * derivative(x[inside])
  }
  values
}

weak_component_rows <- function(test_indices, n_tests, n_components) {
  unlist(lapply(seq_len(n_components), function(component) (component - 1L) * n_tests + test_indices), use.names = FALSE)
}

# Interior tests sit at every grid point whose support fits in the span, for
# each radius. Default radii halve 0.4 * span down to one grid spacing.
build_weak_test_pool <- function(grid, control) {
  span <- diff(range(grid))
  start <- min(grid)
  end <- max(grid)
  radii <- control$test_radii
  if (is.null(radii)) {
    spacing <- span / max(1L, length(grid) - 1L)
    radii <- .4 * span / 2^(0:max(0L, floor(log2(.4 * span / spacing))))
  }
  radii <- sort(unique(radii))
  if (any(radii >= span / 2)) stop("test_radii must be smaller than half the span.")
  # Admit supports that end at the span edges up to roundoff.
  tol <- 32 * .Machine$double.eps * max(1, abs(start), abs(end), span)
  interior <- do.call(rbind, lapply(radii, function(radius) {
    centers <- grid[grid >= start + radius - tol & grid <= end - radius + tol]
    data.frame(center = centers, radius = rep(radius, length(centers)))
  }))
  if (!nrow(interior)) stop("No admissible interior grid centers; refine the working grid.")
  boundary_radii <- if (is.null(control$boundary_radii)) max(radii) else control$boundary_radii
  if (any(boundary_radii >= span / 2)) stop("boundary_radii must be smaller than half the span.")
  boundary <- if (control$include_boundary) {
    do.call(rbind, lapply(boundary_radii, function(radius) {
      offsets <- seq(0, .4 * radius, length.out = control$boundary_count)
      data.frame(center = c(start + offsets, end - offsets), radius = radius)
    }))
  } else {
    data.frame(center = numeric(), radius = numeric())
  }
  list(interior = interior, boundary = boundary, radii = radii, boundary_radii = boundary_radii)
}

evaluate_weak_test_rows <- function(tt, design, order, eta) {
  if (!nrow(design)) {
    return(matrix(0, 0, length(tt)))
  }
  scaled_distance <- (matrix(rep(tt, each = nrow(design)), nrow = nrow(design)) - design$center) / design$radius
  evaluate_bump(scaled_distance, order, eta) / design$radius^order
}

# map %*% V has orthonormal rows spanning the singular directions of the
# row-normalized V above tol times the largest singular value.
orthonormalize_weak_tests <- function(V, tol) {
  if (!nrow(V)) {
    return(matrix(0, 0, 0))
  }
  norms <- sqrt(rowSums(V^2))
  if (any(norms == 0)) stop("Weak supports are unresolved on the initial grid.")
  normalized_tests <- V / norms
  decomposition <- svd(normalized_tests, nu = min(dim(normalized_tests)), nv = 0)
  keep <- which(decomposition$d > max(decomposition$d) * tol)
  sweep(t(decomposition$u[, keep, drop = FALSE]) / decomposition$d[keep], 2, norms, "/")
}

# Orthonormal test combinations on the grid, keeping every numerically
# distinguishable direction (about one per interior grid point). The full span
# is what pins down grid-scale modes of a latent state.
build_test_basis <- function(grid, design, control) {
  list(
    interior = orthonormalize_weak_tests(evaluate_weak_test_rows(grid, design$interior, 0L, control$bump_eta), control$basis_tol),
    boundary = orthonormalize_weak_tests(evaluate_weak_test_rows(grid, design$boundary, 0L, control$bump_eta), control$basis_tol)
  )
}

# ---- ODE model and quadrature residual -------------------------------------
build_symbolic_ode_model <- function(f, n_components, n_parameters) {
  u <- sym_build(lapply(seq_len(n_components), function(j) sym_symbol(paste0("u", j))))
  p <- sym_build(lapply(seq_len(n_parameters), function(j) sym_symbol(paste0("p", j))))
  time <- sym_symbol("t")
  variables <- c(p, u, time)
  rhs <- f(u, p, time)
  if (sym_length(rhs) != n_components) stop("f must return one derivative per observed component.")
  list(
    D = n_components, J = n_parameters,
    rhs = build_fn(rhs, vars = variables),
    jacobian = build_fn(compute_symbolic_jacobian(rhs, c(p, u)), variables)
  )
}

# Golub-Welsch nodes and weights on [-1, 1], mapped to every grid interval.
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

# Stencil width (4 to 16 nodes) follows the fitted GP radius, fixed per grid
# interval so each Gauss cell integrates one polynomial.
choose_lagrange_stencil_size <- function(fit, grid, tt) {
  spacing <- mean(diff(grid))
  n_grid <- length(grid)
  midpoints <- (head(grid, -1L) + tail(grid, -1L)) / 2
  radius <- evaluate_gp_radius((midpoints - fit$origin) / fit$span, fit$radius_coef, fit$radius_bounds)$r * fit$span
  width <- pmin(2L * pmax(2L, pmin(8L, as.integer(round(radius / spacing / 2)))), n_grid)
  width[findInterval(tt, grid, rightmost.closed = TRUE, all.inside = TRUE)]
}

# Centered Lagrange stencil rows. Weights sum to one, so the mean offset cancels.
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

# Per-component linear map from grid values to values at tt. Components with
# data use Lagrange; a latent component uses the GP kriging weights
# K(tt, grid) K(grid, grid)^-1, which suppress its grid-scale modes.
build_state_extensions <- function(fits, factors, grid, tt, latent = integer()) {
  lapply(seq_along(fits), function(component) {
    if (component %in% latent) {
      t(solve_cholesky_system(t(factors[[component]]), t(gp_cross_covariance(fits[[component]], tt, grid))))
    } else {
      build_lagrange_interpolation(grid, tt, choose_lagrange_stencil_size(fits[[component]], grid, tt))
    }
  })
}

# Operator for R = V f(u) + Vp u - B u(endpoints), with the retained tests
# given by maps applied to the test pool.
build_gauss_weak_operator <- function(grid, design, fits, factors, control, maps, latent = integer()) {
  quad <- build_gauss_rule(grid, control$quad_order)
  raw <- evaluate_gauss_test_rows(quad$tt, design, 0L, control$bump_eta)
  n_raw_interior <- nrow(design$interior)
  n_raw_boundary <- nrow(design$boundary)
  n_interior <- nrow(maps$interior)
  n_boundary <- nrow(maps$boundary)
  transform <- function(x) {
    rbind(
      maps$interior %*% x[seq_len(n_raw_interior), , drop = FALSE],
      maps$boundary %*% x[n_raw_interior + seq_len(n_raw_boundary), , drop = FALSE]
    )
  }
  phi <- transform(raw)
  prime <- transform(evaluate_gauss_test_rows(quad$tt, design, 1L, control$bump_eta))
  endpoints <- transform(evaluate_gauss_test_rows(range(grid), design, 0L, control$bump_eta))
  state_extensions <- build_state_extensions(fits, factors, grid, quad$tt, latent)
  boundary_extensions <- build_state_extensions(fits, factors, grid, range(grid), latent)
  weak <- list(
    quad_tt = quad$tt, ni = n_interior, nb = n_boundary, K = n_interior + n_boundary, design = design,
    V = sweep(phi, 2, quad$weights, "*"), Vp = sweep(prime, 2, quad$weights, "*"),
    B = cbind(-endpoints[, 1], endpoints[, 2]),
    gram = tcrossprod(sweep(phi, 2, sqrt(quad$weights), "*")),
    E = state_extensions, Eb = boundary_extensions, mean = vapply(fits, `[[`, numeric(1), "mean")
  )
  # Constant Jacobian blocks, computed once.
  weak$VpE <- lapply(state_extensions, function(E) weak$Vp %*% E)
  weak$BEb <- lapply(boundary_extensions, function(Eb) weak$B %*% Eb)
  weak
}

# Weak residual r (tests x components, component-major) with Jacobians in
# physical parameters (Jp) and grid states (Ju).
evaluate_weak_residual <- function(U, p, weak, model, jacobian = TRUE) {
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
  rhs <- model$rhs(input)
  residual <- weak$V %*% rhs + weak$Vp %*% quadrature_state - weak$B %*% boundary_state
  parameter_jacobian <- NULL
  state_jacobian <- NULL
  if (jacobian) {
    rhs_jacobian <- array(model$jacobian(input), c(n_quadrature, n_components, n_parameters + n_components))
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
        # Scales row q of E by df_a/du_b at quadrature node q.
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

# ---- ODE penalty weights ---------------------------------------------------
# Whitener W with W C W' = I on the numerically nonsingular subspace of C.
# Rank is decided on the correlation matrix, so row rescaling does not matter.
build_covariance_whitener <- function(covariance, tol) {
  if (!nrow(covariance)) {
    return(matrix(0, 0, 0))
  }
  covariance <- (covariance + t(covariance)) / 2
  sd <- sqrt(pmax(diag(covariance), 0))
  active <- which(sd > 0)
  if (!length(active)) stop("Residual covariance has zero rank.")
  correlation <- covariance[active, active, drop = FALSE] / outer(sd[active], sd[active])
  decomposition <- eigen(correlation, symmetric = TRUE)
  if (min(decomposition$values) < -1e-6 * max(decomposition$values)) stop("Residual covariance is not positive semidefinite.")
  keep <- which(decomposition$values > max(decomposition$values) * tol)
  W <- matrix(0, length(keep), nrow(covariance))
  W[, active] <- sweep(
    t(decomposition$vectors[, keep, drop = FALSE]) / sqrt(decomposition$values[keep]),
    2, sd[active], "/"
  )
  W
}

# Per-component ODE precision T / rms^2, which makes the penalty dimensionless.
calculate_ode_scaling <- function(Y, tt, control) {
  rms <- sqrt(colMeans(Y^2))
  span <- diff(range(tt))
  if (!is.null(control$ode_component_scale)) {
    if (!length(control$ode_component_scale) %in% c(1L, ncol(Y))) {
      stop("ode_component_scale must be scalar or have one value per component.")
    }
    rms <- rep_len(control$ode_component_scale, ncol(Y))
  }
  if (!is.finite(span) || span <= 0 || any(!is.finite(rms) | rms <= 0)) {
    stop("RMS/span scaling requires a positive span and positive finite component RMS.")
  }
  precision <- span / rms^2
  if (any(!is.finite(precision) | precision <= 0)) stop("Nonfinite RMS/span precision; rescale input units.")
  list(rms = rms, precision = precision)
}

# Test Gram matrix over all K * D residual rows in component-major order. rows
# selects a principal block (the interior rows during initialization).
build_ode_metric <- function(weak, n_components, rows = seq_len(weak$K * n_components)) {
  kronecker(diag(n_components), weak$gram)[rows, rows, drop = FALSE]
}

# Whitening weights for the ODE residual: interior rows are whitened by the
# metric, boundary rows conditionally on the interior rows. Columns are then
# scaled by the component precision.
build_ode_weights <- function(metric, weak, n_components, control, precision) {
  interior_rows <- weak_component_rows(seq_len(weak$ni), weak$K, n_components)
  boundary_rows <- setdiff(seq_len(weak$K * n_components), interior_rows)
  interior <- build_covariance_whitener(metric[interior_rows, interior_rows, drop = FALSE], control$covariance_tol)
  Wi <- matrix(0, nrow(interior), nrow(metric))
  Wi[, interior_rows] <- interior
  Wb <- matrix(0, 0, nrow(metric))
  if (length(boundary_rows)) {
    whitened_cross <- interior %*% metric[interior_rows, boundary_rows, drop = FALSE]
    conditional <- metric[boundary_rows, boundary_rows, drop = FALSE] - crossprod(whitened_cross)
    # Clip negative roundoff on the Schur complement diagonal.
    diag(conditional) <- pmax(diag(conditional), 0)
    conditional_residual <- matrix(0, length(boundary_rows), nrow(metric))
    conditional_residual[, boundary_rows] <- diag(length(boundary_rows))
    conditional_residual[, interior_rows] <- -crossprod(whitened_cross, interior)
    Wb <- build_covariance_whitener(conditional, control$covariance_tol) %*% conditional_residual
  }
  factors <- rep(sqrt(precision), each = weak$K)
  Wi <- sweep(Wi, 2, factors, "*")
  Wb <- sweep(Wb, 2, factors, "*")
  list(W = rbind(Wi, Wb), Wi = Wi, Wb = Wb)
}
