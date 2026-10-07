# The data, GP prior, and weak ODE terms share one residual and Jacobian.
# Coordinates only change how parameters and states are represented.
make_joint_gp_objective <- function(Y, H, weak, model, weights, noise, lambda,
                                    coordinates, prior_state = NULL,
                                    prior_components = integer(),
                                    observed = seq_len(model$D),
                                    prior_whitened = FALSE) {
  n_grid <- nrow(coordinates$mu)
  n_components <- model$D
  n_parameters <- model$J
  n_observations <- nrow(Y)
  state_columns <- lapply(seq_len(n_components), function(component) {
    n_parameters + (match(component, coordinates$order) - 1L) * n_grid + seq_len(n_grid)
  })
  transform_state_jacobian <- function(block, component) {
    coordinate_map <- coordinates$L[[component]]
    if (length(coordinate_map) == 1L) block * coordinate_map else block %*% coordinate_map
  }

  # These rows are linear and do not depend on the optimization point.
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
    block <- if (prior_whitened) {
      diag(n_grid)
    } else {
      forwardsolve(prior_state$L[[component]], transform_state_jacobian(diag(n_grid), component))
    }
    linear_jacobian[rows, state_columns[[component]]] <- block
  }

  function(theta, jacobian = TRUE) {
    point <- unpack_joint_gp_coordinates(theta, coordinates)
    U <- point$U
    p <- point$p
    data_residual <- sweep(Y - H %*% U[, observed, drop = FALSE], 2, noise, "/")
    prior_residual <- matrix(0, n_grid, length(prior_components))
    for (k in seq_along(prior_components)) {
      component <- prior_components[k]
      prior_residual[, k] <- if (prior_whitened) {
        theta[state_columns[[component]]]
      } else {
        forwardsolve(prior_state$L[[component]], U[, component] - prior_state$mu[, component])
      }
    }
    weak_value <- evaluate_weak_residual(U, p, weak, model, jacobian)
    ode_residual <- as.vector(sqrt(lambda) * weights$W %*% weak_value$r)
    residual <- c(data_residual, prior_residual, ode_residual)
    interior <- seq_len(weights$interior$rank)
    boundary <- seq_len(weights$boundary$rank) + length(interior)
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
    scale <- pmax(sqrt(colSums(current$J^2)), 1e-8)
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
      small_step <- sqrt(sum(step^2)) <= control$xtol * (control$xtol + sqrt(sum(x^2)))
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
# uses unweighted interior residuals; later passes use the chosen ODE metric.
initialize_ode_parameters <- function(U, p0, weak, model, Sigma, control, lower, upper) {
  rows <- weak_component_rows(seq_len(weak$ni), weak$K, model$D)
  p <- p0
  for (iteration in 0:control$init_gls) {
    whitener <- diag(length(rows))
    if (iteration > 0 && control$ode_weighting != "identity") {
      weak_value <- evaluate_weak_residual(U, p, weak, model)
      metric <- build_ode_metric(weak, model$D, control, weak_value$Ju, Sigma, rows)
      whitener <- build_covariance_whitener(metric, control$covariance_tol)$W
    }
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

pack_initial_gp_state <- function(state) {
  coordinates <- state$coordinates
  coordinates$pscale <- numeric()
  pack_joint_gp_coordinates(numeric(), state$mean, coordinates)
}

# Observed-coordinate preparation for the shared objective.
make_observed_gp_objective <- function(Y, H, state, weak, model, weights, noise, lambda) {
  use_prior <- !isFALSE(state$include_gp_prior)
  coordinates <- state$coordinates
  coordinates$pscale <- rep(1, model$J)
  make_joint_gp_objective(Y, H, weak, model, weights, noise, lambda, coordinates, state,
    if (use_prior) seq_len(model$D) else integer(),
    prior_whitened = state$prior_whitened
  )
}
