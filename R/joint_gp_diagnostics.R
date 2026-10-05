# Final diagnostics use the shared residual engine in physical (p,U)
# coordinates, independently of each optimization preconditioner.
diagnose_joint_gp_fit <- function(U, p, Y, H, noise, weak, model, weights, control,
                                  prior_state = NULL, prior_components = integer(),
                                  observed = seq_len(model$D),
                                  lower = rep(-Inf, length(p)),
                                  upper = rep(Inf, length(p)),
                                  tolerance = control$gtol, rank = TRUE,
                                  state_bounds = attr(control, "state_bounds")) {
  n_grid <- nrow(U)
  n_components <- ncol(U)
  n_parameters <- length(p)
  evaluate <- make_joint_gp_objective(Y, H, weak, model, weights, noise, control$lambda,
    coordinates = list(
      mu = matrix(0, n_grid, n_components), L = rep(list(1), n_components),
      order = seq_len(n_components), pscale = rep(1, n_parameters)
    ),
    prior_state = prior_state, prior_components = prior_components, observed = observed
  )
  theta <- c(p, as.vector(U))
  final <- evaluate(theta)
  residual <- final$r
  jacobian <- final$J
  gradient <- final$gradient
  if (is.null(state_bounds)) state_bounds <- validate_joint_gp_bounds(NULL, NULL, n_components, "state")
  lower_bounds <- c(lower, rep(state_bounds$lower, each = n_grid))
  upper_bounds <- c(upper, rep(state_bounds$upper, each = n_grid))
  active_tolerance <- 32 * .Machine$double.eps * pmax(
    1, abs(theta),
    ifelse(is.finite(lower_bounds), abs(lower_bounds), 0), ifelse(is.finite(upper_bounds), abs(upper_bounds), 0)
  )
  projected_gradient <- gradient
  at_lower <- theta <= lower_bounds + active_tolerance & gradient > 0
  at_upper <- theta >= upper_bounds - active_tolerance & gradient < 0
  projected_gradient[at_lower | at_upper] <- 0
  column_norms <- sqrt(colSums(jacobian^2))
  scaled_gradient <- max(abs(projected_gradient) / pmax(column_norms, 1e-8))
  finite <- all(is.finite(c(residual, gradient, column_norms))) && is.finite(sum(residual^2))
  local_rank <- NULL
  if (rank && finite) {
    singular_values <- svd(sweep(jacobian, 2, pmax(column_norms, 1e-8), "/"), nu = 0, nv = 0)$d
    cutoff <- sqrt(.Machine$double.eps)
    retained <- if (length(singular_values) && max(singular_values) > 0) {
      sum(singular_values > max(singular_values) * cutoff)
    } else {
      0L
    }
    local_rank <- list(
      rank = retained, variables = ncol(jacobian), nullity = ncol(jacobian) - retained,
      singular_values = singular_values, relative_tolerance = cutoff,
      interpretation = "local joint residual Jacobian rank; not global identifiability"
    )
  }
  feasible <- all(theta >= lower_bounds - active_tolerance & theta <= upper_bounds + active_tolerance)
  list(
    stationary = finite && feasible && scaled_gradient <= tolerance, feasible = feasible, scaled_gradient = scaled_gradient,
    gradient = gradient, projected_gradient = projected_gradient, column_norms = column_norms,
    tolerance = tolerance, coordinates = "physical", local_rank = local_rank,
    criterion = "Jacobian-column-scaled projected gradient"
  )
}

measure_gp_interpolation_error <- function(fit, tt, obs, H) {
  n_observations <- length(obs)
  if (all(rowSums(H != 0) == 1L)) {
    return(0)
  }
  prediction <- predict_component_gp(fit, c(obs, tt))
  interpolation_difference <- cbind(diag(n_observations), -H)
  bias <- as.vector(interpolation_difference %*% prediction$mean)
  variance <- rowSums((interpolation_difference %*% prediction$Sigma) * interpolation_difference)
  max(sqrt(pmax(0, bias^2 + variance))) / sqrt(fit$noise2)
}

# Report RMS differences in the same metric as the fitted weak residual.
measure_weak_residual_error <- function(change, weights) {
  rms <- function(W) if (nrow(W)) sqrt(mean((W %*% change)^2)) else 0
  c(interior = rms(weights$Wi), boundary = rms(weights$Wb))
}

# Compare the initial weak residual on nested trapezoidal grids.
check_initial_grid_refinement <- function(state, p, weak, model, fits, weights, control) {
  if (control$lambda == 0) {
    return(c(interior = 0, boundary = 0))
  }
  fine_tt <- seq(min(weak$tt), max(weak$tt), length.out = 2 * length(weak$tt) - 1L)
  fine_mean <- do.call(cbind, lapply(fits, function(f) predict_component_gp(f, fine_tt)$mean))
  fine_weak <- build_grid_weak_operator(fine_tt, weak$design, control, weak$maps)
  coarse <- evaluate_weak_residual(state$mean, p, weak, model, FALSE)$r
  fine <- evaluate_weak_residual(fine_mean, p, fine_weak, model, FALSE)$r
  measure_weak_residual_error(coarse - fine, weights)
}

check_joint_gp_accuracy <- function(U, p, state, weak, model, fits, weights, H, control) {
  if (identical(weak$integration, "gp_gauss")) {
    return(check_gauss_quadrature_accuracy(U, p, state, weak, model, fits, weights, H, control))
  }
  fine_tt <- seq(min(weak$tt), max(weak$tt), length.out = 2 * length(weak$tt) - 1L)
  fine_U <- matrix(0, length(fine_tt), ncol(U))
  observation_error <- numeric(ncol(U))
  observations_on_grid <- all(rowSums(H != 0) == 1L)
  for (component in seq_len(ncol(U))) {
    if (isFALSE(state$include_gp_prior)) {
      # A free grid state has the piecewise-linear representation used by H,
      # not a GP-conditioned extension that would reintroduce GP assumptions.
      fine_U[, component] <- build_observation_matrix(weak$tt, fine_tt) %*% U[, component]
      next
    }
    # GP conditional extension of the optimized grid state. Numerical jitter is
    # only on existing grid variables, not a new white-noise signal between them.
    alpha <- solve_cholesky_system(t(state$L[[component]]), U[, component] - state$mu[, component])
    fine_U[, component] <- fits[[component]]$mean + gp_cross_covariance(fits[[component]], fine_tt, weak$tt) %*% alpha
    if (!observations_on_grid) {
      observed <- fits[[component]]$mean + gp_cross_covariance(fits[[component]], fits[[component]]$tt, weak$tt) %*% alpha
      observation_error[component] <- max(abs(observed - H %*% U[, component])) / sqrt(fits[[component]]$noise2)
    }
  }
  fine_weak <- build_grid_weak_operator(fine_tt, weak$design, control, weak$maps)
  # The nested grid retains the fitted values at its original nodes.
  fine_U[seq(1L, length(fine_tt), by = 2L), ] <- U
  coarse <- evaluate_weak_residual(U, p, weak, model, FALSE)$r
  fine <- evaluate_weak_residual(fine_U, p, fine_weak, model, FALSE)$r
  errors <- if (control$lambda > 0) {
    measure_weak_residual_error(coarse - fine, weights)
  } else {
    c(interior = 0, boundary = 0)
  }
  list(
    interpolation = observation_error, weak = errors,
    weak_passed = max(errors) <= control$weak_grid_tol,
    passed = max(observation_error) <= control$grid_tol
  )
}


check_gauss_quadrature_accuracy <- function(U, p, state, weak, model, fits, weights, H, control,
                                            observed = seq_len(model$D),
                                            observation_times = fits[[observed[1L]]]$tt) {
  observation_error <- setNames(numeric(length(observed)), observed)
  observations_on_grid <- all(rowSums(H != 0) == 1L)
  # Use the representation stored with the objective, not a caller's default.
  control$weak_extension <- weak$extension_kind
  if (!observations_on_grid) {
    for (k in seq_along(observed)) {
      component <- observed[k]
      component_control <- control
      component_control$weak_extension <- rep_len(weak$extension_kind, model$D)[component]
      observation_extension <- build_state_extensions(
        fits[component], list(L = state$L[component]), weak$tt, observation_times, component_control
      )[[1]]
      extended_state <- weak$mean[component] + observation_extension %*% (U[, component] - weak$mean[component])
      observation_error[k] <- max(abs(extended_state - H %*% U[, component])) / sqrt(fits[[component]]$noise2)
    }
  }
  fine <- build_gauss_weak_operator(
    weak$tt, weak$design, fits, state,
    control, weak$maps, 2L * weak$quad_order
  )
  coarse_value <- evaluate_gauss_weak_residual(U, p, weak, model)
  fine_value <- evaluate_gauss_weak_residual(U, p, fine, model)
  unweighted_error <- measure_weak_residual_error(coarse_value$r - fine_value$r, weights)
  weak_error <- sqrt(control$lambda) * unweighted_error
  # Compare derivatives in declared, fixed physical-coordinate scales. This
  # does not depend on the optimization coordinates.
  scales <- c(weak$parameter_scale, rep(weak$state_scale, each = nrow(U)))
  coarse_jacobian <- sweep(weights$W %*% cbind(coarse_value$Jp, coarse_value$Ju), 2, scales, "*")
  fine_jacobian <- sweep(weights$W %*% cbind(fine_value$Jp, fine_value$Ju), 2, scales, "*")
  jacobian_relative <- norm(coarse_jacobian - fine_jacobian, "F") / max(norm(fine_jacobian, "F"), .Machine$double.eps)
  gram_whitener <- build_covariance_whitener(weak$gram, control$covariance_tol)$W
  gram_error <- norm(gram_whitener %*% (fine$gram - weak$gram) %*% t(gram_whitener), "2")
  quadrature_passed <- all(is.finite(c(weak_error, jacobian_relative, gram_error))) &&
    max(weak_error) <= control$weak_grid_tol &&
    jacobian_relative <= control$weak_design_quad_tol && gram_error <= control$weak_quad_gram_tol
  list(
    interpolation = observation_error,
    observed = observed, weak = weak_error, weak_unweighted = unweighted_error,
    jacobian_relative = jacobian_relative, gram_error = gram_error,
    quadrature_passed = quadrature_passed,
    weak_passed = quadrature_passed,
    passed = all(is.finite(observation_error)) && max(observation_error) <= control$grid_tol,
    quad_order = weak$quad_order, check_order = fine$quad_order,
    check_points = length(fine$quad_tt),
    quadrature_points = length(weak$quad_tt), state_variables = length(U)
  )
}

format_quadrature_warning <- function(check, control) {
  values <- c(residual = max(check$weak), Jacobian = check$jacobian_relative, Gram = check$gram_error)
  tolerances <- c(control$weak_grid_tol, control$weak_design_quad_tol, control$weak_quad_gram_tol)
  failed <- !is.finite(values) | values > tolerances
  details <- sprintf(
    "%s discrepancy %.3g (tolerance %.3g)", names(values)[failed],
    values[failed], tolerances[failed]
  )
  paste0(
    "Final weak quadrature check failed: ", paste(details, collapse = "; "),
    ". Inspect fit$diagnostics$final_grid.",
    if (control$weak_quad_order < 32L) {
      " A higher control$weak_quad_order refines integration without adding state variables."
    } else {
      ""
    }
  )
}
