# Weak residuals, symbolic derivatives, and residual weighting.

build_symbolic_ode_model <- function(f, n_components, n_parameters, em_order) {
  u <- sym_build(lapply(seq_len(n_components), function(j) sym_symbol(paste0("u", j))))
  p <- sym_build(lapply(seq_len(n_parameters), function(j) sym_symbol(paste0("p", j))))
  time <- sym_symbol("t")
  variables <- c(p, u, time)
  rhs <- f(u, p, time)
  if (sym_length(rhs) != n_components) stop("f must return one derivative per observed component.")
  jets <- list(rhs)
  n_derivatives <- if (em_order == 4L) 4L else if (em_order == 2L) 2L else 1L
  if (n_derivatives > 1L) {
    for (j in 2:n_derivatives) {
      jets[[j]] <- compute_symbolic_total_time_deriv(jets[[j - 1L]], u, rhs, time)
    }
  }
  parameters_and_state <- c(p, u)
  model <- list(
    D = n_components, J = n_parameters,
    jet = lapply(jets, build_fn, vars = variables),
    jet_jac = lapply(jets, function(z) build_fn(compute_symbolic_jacobian(z, parameters_and_state), variables))
  )
  model
}

# Physical trapezoidal weights, including irregular observation grids. The
# diagnostic grid and observation grid must use the SAME integral convention.
trapezoidal_weights <- function(tt) {
  interval_widths <- diff(tt)
  c(interval_widths[1], head(interval_widths, -1L) + tail(interval_widths, -1L), interval_widths[length(interval_widths)]) / 2
}

# E ||A U||^2 / nrow(A), conditional on fitted GP hyperparameters. B maps
# ORIGINAL measurement errors, not independent pseudo-observations on a dense
# grid. No numerical GP jitter is added to either stochastic quantity here.
posterior_fourier_moments <- function(A, mean, Sigma, B, noise2) {
  mean_square <- mean(as.vector(A %*% mean)^2)
  terms <- (A %*% Sigma) * A
  posterior_variance <- sum(terms) / nrow(A)
  roundoff <- 1e-8 * max(abs(diag(Sigma))) * sum(A^2) / nrow(A)
  if (posterior_variance < -max(roundoff, .Machine$double.eps)) {
    stop("Posterior radius diagnostic encountered a negative covariance contribution.")
  }
  posterior_variance <- max(0, posterior_variance)
  noise_variance <- noise2 * sum(B^2) / nrow(B)
  total <- mean_square + posterior_variance
  list(
    mean_square = mean_square, posterior_variance = posterior_variance,
    expected_square = total, noise_variance = noise_variance,
    rms_ratio = if (noise_variance > 0) sqrt(total / noise_variance) else Inf
  )
}

place_weak_test_centers <- function(tt, radius, control, dense = FALSE) {
  span <- diff(range(tt))
  count <- max(3L, min(
    control$weak_count,
    ceiling((span - 2 * radius) / (radius / if (dense) 4 else 2)) + 1L
  ))
  seq(min(tt) + radius, max(tt) - radius, length.out = count)
}

# Posterior version of the high-frequency proxy used by
# find_min_radius_int_error. Both Fourier phases are retained to avoid a zero
# caused merely by phase alignment. Its carrier frequency is tied to the
# ORIGINAL observation count, not the number of interpolated grid values.
build_radius_test_operators <- function(tt, diagnostic_tt, radius, control) {
  centers <- place_weak_test_centers(tt, radius, control, dense = TRUE)
  design <- data.frame(center = centers, radius = radius)
  grid_weights <- trapezoidal_weights(diagnostic_tt)
  observation_weights <- trapezoidal_weights(tt)
  raw <- evaluate_weak_test_rows(diagnostic_tt, design, 0L, control$bump_eta)
  norms <- sqrt(as.vector(raw^2 %*% grid_weights))
  if (any(!is.finite(norms) | norms <= 0)) stop("Unresolved radius diagnostic supports.")
  wavenumber <- max(1L, floor(length(tt) / 3L) - 1L)
  phase <- function(t) 2 * pi * wavenumber * (t - min(tt)) / diff(range(tt))
  operator <- function(values, times, weights) {
    values <- values / norms
    rbind(
      sweep(values, 2, weights * sin(phase(times)), "*"),
      sweep(values, 2, weights * cos(phase(times)), "*")
    )
  }
  list(
    A = operator(raw, diagnostic_tt, grid_weights),
    B = operator(evaluate_weak_test_rows(tt, design, 0L, control$bump_eta), tt, observation_weights),
    centers = centers, wavenumber = wavenumber
  )
}

select_posterior_test_radius <- function(tt, fits, control) {
  # A dedicated diagnostic grid makes physical test supports invariant to
  # working-grid choice. Its resolution can be checked independently.
  size <- max(length(tt), control$weak_radius_grid)
  diagnostic_tt <- seq(min(tt), max(tt), length.out = size)
  span <- diff(range(tt))
  smallest <- 4 * span / (size - 1L)
  largest <- 0.4 * span
  smallest <- min(smallest, largest)
  steps <- floor(log(largest / smallest) / log(control$weak_radius_ratio))
  radii <- sort(unique(c(smallest * control$weak_radius_ratio^(0:steps), largest)))
  predictions <- lapply(fits, predict_component_gp, tt = diagnostic_tt)
  operators <- lapply(radii, function(r) build_radius_test_operators(tt, diagnostic_tt, r, control))
  records <- list()
  index <- 0L
  for (component in seq_along(fits)) {
    for (j in seq_along(radii)) {
      prediction <- predictions[[component]]
      operator <- operators[[j]]
      moments <- posterior_fourier_moments(operator$A, prediction$mean, prediction$Sigma, operator$B, fits[[component]]$noise2)
      index <- index + 1L
      records[[index]] <- data.frame(
        component = component, radius = radii[j],
        mean_square = moments$mean_square, posterior_variance = moments$posterior_variance,
        expected_square = moments$expected_square, noise_variance = moments$noise_variance,
        rms_ratio = moments$rms_ratio, centers = length(operator$centers),
        passed = is.finite(moments$rms_ratio) && moments$rms_ratio <= control$weak_radius_tolerance
      )
    }
  }
  profile <- do.call(rbind, records)
  per_component <- vapply(seq_along(fits), function(component) {
    rows <- profile[profile$component == component, ]
    # An isolated Fourier notch must not admit an otherwise unresolved family.
    passing_rows <- which(as.logical(rev(cumprod(rev(rows$passed)))))
    if (length(passing_rows)) rows$radius[passing_rows[1L]] else NA_real_
  }, numeric(1))
  selection <- list(
    method = "posterior", profile = profile, candidates = radii,
    component_minimum = per_component, diagnostic_tt = diagnostic_tt,
    observation_points = length(tt), wavenumber = operators[[1L]]$wavenumber,
    tolerance = control$weak_radius_tolerance, includes_posterior_covariance = TRUE,
    noise_on_observation_grid = TRUE, hyperparameters_conditioned_on = TRUE,
    criterion = "posterior Fourier-proxy RMS / original-observation noise RMS",
    passed = all(is.finite(per_component))
  )
  if (!selection$passed) {
    failed_components <- which(!is.finite(per_component))
    condition <- structure(
      list(message = paste0(
        "No posterior/noise-floor weak radius is admissible for component(s) ",
        paste(failed_components, collapse = ", "), ". Inspect the radius profile, GP/noise fit and observation resolution; ",
        "no radius was silently substituted."
      ), call = NULL, radius_selection = selection),
      class = c("wendygp_radius_error", "error", "condition")
    )
    stop(condition)
  }
  # The current residual shares one test bank across components, so require the
  # entire retained family to pass for EVERY component, not only a pooled average.
  selection$minimum <- max(per_component)
  selection$radii <- radii[radii >= selection$minimum]
  selection
}

build_weak_test_design <- function(tt, fits, control) {
  span <- diff(range(tt))
  spacing <- min(diff(tt))
  radii <- control$weak_radii
  automatic <- is.null(radii)
  selection <- list(method = if (automatic) control$weak_radius_method else "explicit")
  if (automatic && control$weak_radius_method == "posterior") {
    selection <- select_posterior_test_radius(tt, fits, control)
    radii <- selection$radii
  } else if (automatic) {
    gp_radii <- vapply(fits, function(g) {
      g$span * stats::median(
        evaluate_gp_radius(seq(0, 1, length.out = 25), g$radius_coef, g$radius_bounds)$r
      )
    }, numeric(1))
    radii <- as.vector(outer(gp_radii, control$weak_factors))
    # Physical supports do not change when the working grid is refined.
    radii <- pmin(0.4 * span, pmax(4 * spacing, radii))
    # Snap proposals to a shared geometric family. Nearly identical radii from
    # different components must not create a forest of derivative-like BL modes.
    family <- sort(unique(pmin(0.4 * span, pmax(4 * spacing, span * 2^(-6:-1)))))
    radii <- sort(unique(vapply(radii, function(r) family[which.min(abs(log(family / r)))], numeric(1))))
  }
  radii <- sort(unique(radii))
  if (any(radii >= span / 2)) stop("weak_radii must be smaller than half the time span.")
  start <- min(tt)
  end <- max(tt)
  interior <- do.call(rbind, lapply(radii, function(r) {
    data.frame(center = place_weak_test_centers(tt, r, control,
      dense = automatic && control$weak_radius_method == "posterior"
    ), radius = r)
  }))
  boundary_radii <- control$bl_radii
  if (is.null(boundary_radii)) {
    boundary_radii <-
      if (automatic && control$weak_radius_method == "posterior") max(radii) else radii
  }
  if (any(boundary_radii >= span / 2)) stop("bl_radii must be smaller than half the time span.")
  boundary_radii <- sort(unique(boundary_radii))
  boundary <- if (control$include_bl) {
    do.call(rbind, lapply(boundary_radii, function(r) {
      offset <- seq(0, 0.4 * r, length.out = control$bl_count)
      data.frame(center = c(start + offset, end - offset), radius = r)
    }))
  } else {
    data.frame(center = numeric(), radius = numeric())
  }
  list(
    interior = interior, boundary = boundary, radii = radii,
    radius_selection = selection, boundary_radii = boundary_radii
  )
}

evaluate_weak_test_rows <- function(tt, design, order, eta) {
  if (!nrow(design)) {
    return(matrix(0, 0, length(tt)))
  }
  scaled_distance <- (matrix(rep(tt, each = nrow(design)), nrow = nrow(design)) - design$center) / design$radius
  evaluate_bump(scaled_distance, order, eta) / design$radius^order
}

orthonormalize_weak_tests <- function(V, tol, spectrum = FALSE) {
  if (!nrow(V)) {
    return(matrix(0, 0, 0))
  }
  norms <- sqrt(rowSums(V^2))
  if (any(norms == 0)) stop("Weak supports are unresolved on the initial grid.")
  normalized_tests <- V / norms
  decomposition <- svd(normalized_tests, nu = min(dim(normalized_tests)), nv = 0)
  keep <- which(decomposition$d > max(decomposition$d) * tol)
  map <- sweep(t(decomposition$u[, keep, drop = FALSE]) / decomposition$d[keep], 2, norms, "/")
  if (spectrum) {
    list(
      map = map, singular_values = decomposition$d,
      retained = keep
    )
  } else {
    map
  }
}

build_grid_weak_operator <- function(tt, design, control, maps = NULL) {
  n_grid <- length(tt)
  spacing <- mean(diff(tt))
  weights <- rep(spacing, n_grid)
  weights[c(1, n_grid)] <- spacing / 2
  interior_tests <- evaluate_weak_test_rows(tt, design$interior, 0L, control$bump_eta)
  boundary_tests <- evaluate_weak_test_rows(tt, design$boundary, 0L, control$bump_eta)
  if (is.null(maps)) {
    maps <- list(
      interior = orthonormalize_weak_tests(interior_tests, control$basis_tol),
      boundary = orthonormalize_weak_tests(boundary_tests, control$basis_tol)
    )
  }
  V <- rbind(maps$interior %*% interior_tests, maps$boundary %*% boundary_tests)
  Vp <- rbind(
    maps$interior %*% evaluate_weak_test_rows(tt, design$interior, 1L, control$bump_eta),
    maps$boundary %*% evaluate_weak_test_rows(tt, design$boundary, 1L, control$bump_eta)
  )
  n_interior <- nrow(maps$interior)
  n_boundary <- nrow(maps$boundary)
  endpoints <- lapply(c(tt[1], tt[n_grid]), function(t) {
    raw <- vapply(0:4, function(j) {
      as.vector(evaluate_weak_test_rows(
        t, design$boundary, j, control$bump_eta
      ))
    }, numeric(nrow(design$boundary)))
    if (!n_boundary) matrix(0, 0, 5) else maps$boundary %*% raw
  })
  V <- sweep(V, 2, weights, "*")
  Vp <- sweep(Vp, 2, weights, "*")
  if (n_boundary) {
    Vp[n_interior + seq_len(n_boundary), 1] <- Vp[n_interior + seq_len(n_boundary), 1] + endpoints[[1]][, 1]
    Vp[n_interior + seq_len(n_boundary), n_grid] <- Vp[n_interior + seq_len(n_boundary), n_grid] - endpoints[[2]][, 1]
  }
  list(
    V = V, Vp = Vp, ni = n_interior, nb = n_boundary, K = n_interior + n_boundary, tt = tt,
    endpoints = endpoints, maps = maps, design = design, h = spacing,
    em_order = control$em_order
  )
}

evaluate_weak_residual <- function(U, p, weak, model, jacobian = TRUE) {
  if (identical(weak$integration, "gp_gauss")) {
    return(evaluate_gauss_weak_residual(U, p, weak, model, jacobian))
  }
  n_grid <- nrow(U)
  n_components <- model$D
  n_parameters <- model$J
  n_tests <- weak$K
  input <- rbind(matrix(p, n_parameters, n_grid), t(U), weak$tt)
  rhs <- model$jet[[1]](input)
  residual <- weak$V %*% rhs + weak$Vp %*% U
  parameter_jacobian <- NULL
  state_jacobian <- NULL
  if (jacobian) {
    rhs_jacobian <- array(model$jet_jac[[1]](input), c(n_grid, n_components, n_parameters + n_components))
    rhs_parameter_jacobian <- matrix(
      rhs_jacobian[, , seq_len(n_parameters), drop = FALSE],
      n_grid, n_components * n_parameters
    )
    parameter_jacobian <- matrix(
      weak$V %*% rhs_parameter_jacobian, n_tests * n_components, n_parameters
    )
    state_jacobian <- matrix(0, n_tests * n_components, n_grid * n_components)
    # Block (a, b) differentiates residual component a by state component b.
    for (a in seq_len(n_components)) {
      for (b in seq_len(n_components)) {
        rows <- (a - 1L) * n_tests + seq_len(n_tests)
        columns <- (b - 1L) * n_grid + seq_len(n_grid)
        block <- sweep(weak$V, 2, rhs_jacobian[, a, n_parameters + b], "*")
        if (a == b) block <- block + weak$Vp
        state_jacobian[rows, columns] <- block
      }
    }
  }
  if (weak$nb && weak$em_order) {
    boundary_rows <- weak$ni + seq_len(weak$nb)
    for (side in 1:2) {
      endpoint_index <- if (side == 1L) 1L else n_grid
      sign <- if (side == 1L) -1 else 1
      endpoint_input <- matrix(c(p, U[endpoint_index, ], weak$tt[endpoint_index]), ncol = 1L)
      endpoint_derivatives <- lapply(model$jet, function(fn) as.vector(fn(endpoint_input)))
      correction <- matrix(0, 5L, n_components)
      correction[1:3, ] <- -weak$h^2 / 12 * g_coeffs(endpoint_derivatives, U[endpoint_index, ], 1L)
      if (weak$em_order == 4L) {
        correction <- correction + weak$h^4 / 720 * g_coeffs(endpoint_derivatives, U[endpoint_index, ], 3L)
      }
      residual[boundary_rows, ] <- residual[boundary_rows, ] + sign * weak$endpoints[[side]] %*% correction
      if (jacobian) {
        endpoint_jacobians <- lapply(model$jet_jac, function(fn) as.vector(fn(endpoint_input)))
        state_derivative <- as.vector(cbind(matrix(0, n_components, n_parameters), diag(n_components)))
        correction_jacobian <- matrix(0, 5L, n_components * (n_parameters + n_components))
        correction_jacobian[1:3, ] <- -weak$h^2 / 12 * g_coeffs(endpoint_jacobians, state_derivative, 1L)
        if (weak$em_order == 4L) {
          correction_jacobian <- correction_jacobian + weak$h^4 / 720 * g_coeffs(endpoint_jacobians, state_derivative, 3L)
        }
        boundary_jacobian <- array(
          sign * weak$endpoints[[side]] %*% correction_jacobian,
          c(weak$nb, n_components, n_parameters + n_components)
        )
        for (a in seq_len(n_components)) {
          rows <- (a - 1L) * n_tests + boundary_rows
          parameter_jacobian[rows, ] <- parameter_jacobian[rows, , drop = FALSE] + matrix(
            boundary_jacobian[, a, seq_len(n_parameters), drop = FALSE],
            weak$nb, n_parameters
          )
          for (b in seq_len(n_components)) {
            columns <- (b - 1L) * n_grid + endpoint_index
            state_jacobian[rows, columns] <- state_jacobian[rows, columns] +
              boundary_jacobian[, a, n_parameters + b]
          }
        }
      }
    }
  }
  list(r = as.vector(residual), Jp = parameter_jacobian, Ju = state_jacobian)
}

# Rank decisions are made on the correlation matrix. They are invariant under
# independent nonzero row rescalings, including sign flips. Discarded modes
# define an explicitly reported projected penalty, never a hidden raw ridge.
build_covariance_whitener <- function(covariance, tol) {
  if (!nrow(covariance)) {
    return(list(W = matrix(0, 0, 0), rank = 0L, discarded = 0L))
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
  list(
    W = W, rank = length(keep), discarded = nrow(covariance) - length(keep),
    smallest_retained = min(decomposition$values[keep]), tolerance = tol
  )
}

build_weak_block_weights <- function(covariance, weak, n_components, tol) {
  interior_rows <- unlist(lapply(seq_len(n_components), function(component) (component - 1L) * weak$K + seq_len(weak$ni)))
  boundary_rows <- setdiff(seq_len(nrow(covariance)), interior_rows)
  interior <- build_covariance_whitener(covariance[interior_rows, interior_rows, drop = FALSE], tol)
  Wi <- matrix(0, interior$rank, nrow(covariance))
  Wi[, interior_rows] <- interior$W
  if (length(boundary_rows)) {
    cross_covariance <- covariance[interior_rows, boundary_rows, drop = FALSE]
    whitened_cross_covariance <- interior$W %*% cross_covariance
    conditional <- covariance[boundary_rows, boundary_rows, drop = FALSE] - crossprod(whitened_cross_covariance)
    # Tiny negative roundoff on a Schur diagonal is not a covariance ridge.
    diag(conditional) <- pmax(diag(conditional), 0)
    boundary <- build_covariance_whitener(conditional, tol)
    conditional_residual_map <- matrix(0, length(boundary_rows), nrow(covariance))
    conditional_residual_map[, boundary_rows] <- diag(length(boundary_rows))
    conditional_residual_map[, interior_rows] <- -crossprod(whitened_cross_covariance, interior$W)
    Wb <- boundary$W %*% conditional_residual_map
  } else {
    boundary <- list(rank = 0L, discarded = 0L)
    conditional <- matrix(0, 0, 0)
    Wb <- matrix(0, 0, nrow(covariance))
  }
  list(
    W = rbind(Wi, Wb), Wi = Wi, Wb = Wb, interior = interior, boundary = boundary,
    conditional_covariance = conditional, ii = interior_rows, bb = boundary_rows
  )
}

propagate_state_covariance <- function(Ju, Sigma, n_grid, n_components) {
  covariance <- matrix(0, nrow(Ju), nrow(Ju))
  for (component in seq_len(n_components)) {
    component_jacobian <- Ju[, (component - 1L) * n_grid + seq_len(n_grid), drop = FALSE]
    covariance <- covariance + component_jacobian %*% Sigma[[component]] %*% t(component_jacobian)
  }
  (covariance + t(covariance)) / 2
}

calculate_ode_scaling <- function(Y, tt, control) {
  units <- control$ode_units
  if (control$ode_weighting != "test_gram") units <- "raw"
  rms <- sqrt(colMeans(Y^2))
  span <- diff(range(tt))
  if (!is.null(control$ode_component_scale)) {
    if (!length(control$ode_component_scale) %in% c(1L, ncol(Y))) {
      stop("ode_component_scale must be scalar or have one value per component.")
    }
    rms <- rep_len(control$ode_component_scale, ncol(Y))
  }
  precision <- rep(1, ncol(Y))
  if (units == "rms_span") {
    if (!is.finite(span) || span <= 0 || any(!is.finite(rms) | rms <= 0)) {
      stop("RMS/span scaling requires a positive span and positive finite component RMS.")
    }
    precision <- span / rms^2
    if (any(!is.finite(precision) | precision <= 0)) stop("Nonfinite RMS/span precision; rescale input units.")
  }
  list(
    units = units, rms = rms, span = span, precision = precision,
    scale_source = if (is.null(control$ode_component_scale)) "observed_rms" else "supplied",
    centered = FALSE, frozen = TRUE
  )
}

# Scale the original retained whiteners, never re-truncate a unit-scaled Gram.
scale_weak_weights <- function(weights, precision, n_tests) {
  stopifnot(all(is.finite(precision) & precision > 0), ncol(weights$W) == n_tests * length(precision))
  if (all(precision == 1)) {
    return(weights)
  }
  factors <- rep(sqrt(precision), each = n_tests)
  scaled_weights <- weights
  for (block in c("W", "Wi", "Wb")) {
    scaled_weights[[block]] <- sweep(weights[[block]], 2, factors, "*")
  }
  scaled_weights$interior$W <- sweep(weights$interior$W, 2, factors[weights$ii], "*")
  if (length(weights$bb)) {
    scaled_weights$boundary$W <- sweep(weights$boundary$W, 2, factors[weights$bb], "*")
    scaled_weights$conditional_covariance <- weights$conditional_covariance /
      outer(factors[weights$bb], factors[weights$bb])
  }
  scaled_weights
}

prepare_ode_metric <- function(metric, weak, Y, tt, control) {
  scaling <- calculate_ode_scaling(Y, tt, control)
  weights <- build_ode_weights(metric, weak, ncol(Y), control)
  factors <- rep(sqrt(scaling$precision), each = weak$K)
  list(
    Omega = metric / outer(factors, factors), raw = metric,
    weights = scale_weak_weights(weights, scaling$precision, weak$K), scaling = scaling
  )
}

# The matrix traditionally stored as Omega is a residual METRIC for the two
# non-GP modes, not an estimate of residual sampling covariance. All rows use
# the existing component-major order. Optional rows avoids full propagation
# during interior-only initialization and gives consistent principal blocks.
build_ode_metric <- function(weak, n_components, control, Ju = NULL, Sigma = NULL, rows = NULL) {
  mode <- control$ode_weighting
  if (is.null(rows)) rows <- seq_len(weak$K * n_components)
  if (mode == "gp_delta") {
    return(propagate_state_covariance(Ju[rows, , drop = FALSE], Sigma, length(weak$tt), n_components))
  }
  if (mode == "identity") {
    return(diag(length(rows)))
  }
  if (mode != "test_gram") stop("Unknown ODE weighting mode.")
  if (identical(weak$integration, "gp_gauss")) {
    metric <- kronecker(diag(n_components), weak$gram)
    return(metric[rows, rows, drop = FALSE])
  }
  # V already contains trapezoidal weights: V V' would incorrectly use dt^2.
  # Recover Phi sqrt(w), giving integral Phi_a(t) Phi_b(t) dt, with boundary
  # tests included. Endpoint/EM derivative corrections belong to the residual,
  # not to this test-function geometry.
  quadrature_weights <- rep(weak$h, length(weak$tt))
  quadrature_weights[c(1L, length(quadrature_weights))] <- weak$h / 2
  test_gram <- tcrossprod(sweep(weak$V, 2, sqrt(quadrature_weights), "/"))
  metric <- kronecker(diag(n_components), test_gram)
  metric[rows, rows, drop = FALSE]
}

build_ode_weights <- function(metric, weak, n_components, control) {
  if (control$ode_weighting != "identity") {
    return(build_weak_block_weights(metric, weak, n_components, control$covariance_tol))
  }
  # No eigenanalysis, row rescaling, or rank truncation in ordinary LS.
  interior_rows <- unlist(lapply(seq_len(n_components), function(component) (component - 1L) * weak$K + seq_len(weak$ni)))
  boundary_rows <- setdiff(seq_len(weak$K * n_components), interior_rows)
  identity <- diag(weak$K * n_components)
  Wi <- identity[interior_rows, , drop = FALSE]
  Wb <- identity[boundary_rows, , drop = FALSE]
  rank_info <- function(n) {
    list(
      rank = n, discarded = 0L,
      smallest_retained = if (n) 1 else NA_real_, tolerance = 0
    )
  }
  list(
    W = rbind(Wi, Wb), Wi = Wi, Wb = Wb, interior = rank_info(length(interior_rows)),
    boundary = rank_info(length(boundary_rows)), conditional_covariance = diag(length(boundary_rows)),
    ii = interior_rows, bb = boundary_rows
  )
}
