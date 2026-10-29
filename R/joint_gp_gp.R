# Component GPs: nonstationary Matern 5/2 kernel, marginal-likelihood fit,
# prediction, and the GP prior on the working grid.

factor_gp_covariance <- function(covariance, jitter = 0) {
  covariance <- (covariance + t(covariance)) / 2
  added <- jitter * max(diag(covariance))
  factor <- chol(covariance + diag(added, nrow(covariance)))
  list(R = factor, added = added)
}

solve_cholesky_system <- function(R, b) backsolve(R, forwardsolve(t(R), b))

evaluate_gp_radius <- function(x, coefficients, bounds) {
  # Keep a locations-by-coefficients matrix even for a single location.
  cosine_basis <- matrix(vapply(
    0:(length(coefficients) - 1L), function(j) cos(j * pi * x),
    numeric(length(x))
  ), nrow = length(x), ncol = length(coefficients))
  bounded_values <- stats::plogis(as.vector(cosine_basis %*% coefficients))
  width <- diff(log(bounds))
  list(
    r = exp(log(bounds[1]) + width * bounded_values),
    derivative = cosine_basis * as.vector(width * bounded_values * (1 - bounded_values))
  )
}

# Paciorek-Schervish nonstationary Matern 5/2 correlation.
evaluate_gp_correlation <- function(x, coefficients, bounds, derivatives = FALSE) {
  radius <- evaluate_gp_radius(x, coefficients, bounds)
  r <- radius$r
  r2 <- r^2
  mean_r2 <- outer(r2, r2, "+") / 2
  z <- sqrt(5) * abs(outer(x, x, "-")) / sqrt(mean_r2)
  prefactor <- sqrt(outer(r, r) / mean_r2)
  decay <- exp(-z)
  R <- prefactor * (1 + z + z^2 / 3) * decay
  # radius$derivative is d log(r) / d coefficient. This form never divides by
  # the distance, so repeated times and the diagonal stay finite.
  radial <- if (derivatives) prefactor * decay * z^2 * (1 + z) / 6
  dR <- if (derivatives) {
    lapply(seq_along(coefficients), function(j) {
      dlogr <- radius$derivative[, j]
      dlog_mean_r2 <- outer(r2 * dlogr, r2 * dlogr, "+") / mean_r2
      R * (outer(dlogr, dlogr, "+") - dlog_mean_r2) / 2 + radial * dlog_mean_r2
    })
  } else {
    NULL
  }
  list(R = R, dR = dR)
}

fit_component_gp <- function(tt, y, noise_sd, control) {
  origin <- min(tt)
  span <- diff(range(tt))
  x <- (tt - origin) / span
  n_obs <- length(y)
  center <- mean(y)
  scale <- stats::sd(y)
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps) * max(abs(y))) {
    stop("A component has essentially constant observations; its GP scale is not identifiable.")
  }
  y_scaled <- (y - center) / scale
  bounds <- if (is.null(control$gp_radius_bounds)) {
    c(max(0.02, min(0.2, 2 * stats::median(diff(x)))), 2)
  } else {
    control$gp_radius_bounds / span
  }
  n_coef <- control$gp_nonstationary_terms + 1L
  penalty_weights <- (0:(n_coef - 1L))^4 * control$gp_radius_penalty
  known_noise <- !is.null(noise_sd)
  noise2 <- if (known_noise) (noise_sd / scale)^2 else NULL
  evaluate <- function(par) {
    coef <- par[seq_len(n_coef)]
    variance_scale <- exp(par[n_coef + 1L])
    correlation <- evaluate_gp_correlation(x, coef, bounds, TRUE)
    covariance <- if (known_noise) {
      variance_scale * correlation$R + diag(noise2, n_obs)
    } else {
      correlation$R + diag(variance_scale, n_obs)
    }
    # Reject numerically indefinite trial covariances instead of failing.
    R <- tryCatch(chol(covariance), error = function(error) NULL)
    if (is.null(R)) {
      return(list(value = 1e50, gradient = rep(0, length(par))))
    }
    precision <- solve_cholesky_system(R, diag(n_obs))
    precision_sums <- rowSums(precision)
    mu <- sum(precision_sums * y_scaled) / sum(precision_sums)
    residual <- y_scaled - mu
    alpha <- as.vector(precision %*% residual)
    tau2 <- if (known_noise) variance_scale else max(sum(residual * alpha) / n_obs, .Machine$double.eps)
    data_term <- if (known_noise) {
      sum(residual * alpha) / 2
    } else {
      n_obs * log(tau2) / 2
    }
    value <- sum(log(diag(R))) + data_term
    value <- value + sum(penalty_weights * coef^2) / 2
    score <- (precision - tcrossprod(alpha) / if (known_noise) 1 else tau2) / 2
    gradient <- vapply(correlation$dR, function(dR) sum(score * dR) * if (known_noise) variance_scale else 1, numeric(1))
    gradient <- c(
      gradient + penalty_weights * coef,
      if (known_noise) sum(score * (variance_scale * correlation$R)) else variance_scale * sum(diag(score))
    )
    list(
      value = value, gradient = gradient, mean = mu, tau2 = tau2,
      noise2 = if (known_noise) noise2 else tau2 * variance_scale, R = correlation$R
    )
  }
  starts <- lapply(seq_len(control$gp_restarts), function(j) {
    fraction <- seq(0.3, 0.8, length.out = control$gp_restarts)[j]
    c(stats::qlogis(fraction), rep(0, n_coef - 1L), if (known_noise) 0 else log(0.03))
  })
  fits <- lapply(starts, function(init) {
    # The marginal likelihood is not a sum of squares, so use nlminb, not LM.
    last_parameters <- NULL
    last_value <- NULL
    evaluate_at <- function(par) {
      if (!identical(par, last_parameters)) {
        last_value <<- evaluate(par)
        last_parameters <<- par
      }
      last_value
    }
    stats::nlminb(init,
      objective = function(par) evaluate_at(par)$value,
      gradient = function(par) evaluate_at(par)$gradient,
      lower = c(rep(-8, n_coef), -18), upper = c(rep(8, n_coef), 10),
      control = list(
        iter.max = control$gp_maxit, eval.max = 4 * control$gp_maxit,
        rel.tol = control$ftol, x.tol = control$xtol
      )
    )
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "objective"))]]
  value <- evaluate(best$par)
  observation_covariance <- value$tau2 * scale^2 * value$R + diag(value$noise2 * scale^2, n_obs)
  observation_factor <- factor_gp_covariance(observation_covariance)$R
  structure(list(
    tt = tt, x = x, y = y, origin = origin, span = span,
    mean = center + scale * value$mean, tau2 = scale^2 * value$tau2,
    noise2 = scale^2 * value$noise2, radius_coef = best$par[seq_len(n_coef)],
    radius_bounds = bounds, Ryy = observation_factor,
    alpha = solve_cholesky_system(observation_factor, y - center - scale * value$mean),
    convergence = best$convergence, message = best$message, objective = best$objective
  ), class = "wendygp_gp")
}

# Posterior mean and prior covariance K at tt.
predict_component_gp <- function(fit, tt) {
  x <- (tt - fit$origin) / fit$span
  all_locations <- c(fit$x, x)
  n_obs <- length(fit$x)
  prediction_indices <- n_obs + seq_along(x)
  R <- evaluate_gp_correlation(all_locations, fit$radius_coef, fit$radius_bounds)$R
  cross_covariance <- fit$tau2 * R[prediction_indices, seq_len(n_obs), drop = FALSE]
  list(
    mean = as.vector(fit$mean + cross_covariance %*% fit$alpha),
    K = fit$tau2 * R[prediction_indices, prediction_indices, drop = FALSE]
  )
}

gp_cross_covariance <- function(fit, t1, t2) {
  x <- (t1 - fit$origin) / fit$span
  y <- (t2 - fit$origin) / fit$span
  rx <- evaluate_gp_radius(x, fit$radius_coef, fit$radius_bounds)$r
  ry <- evaluate_gp_radius(y, fit$radius_coef, fit$radius_bounds)$r
  mean_r2 <- outer(rx^2, ry^2, "+") / 2
  z <- sqrt(5) * abs(outer(x, y, "-")) / sqrt(mean_r2)
  fit$tau2 * sqrt(outer(rx, ry) / mean_r2) * (1 + z + z^2 / 3) * exp(-z)
}

# GP prior (mu, L) and posterior mean on the grid. Unbounded components are
# optimized in prior coordinates U = mu + L z.
prepare_gp_state <- function(fits, tt, control) {
  predictions <- lapply(fits, predict_component_gp, tt = tt)
  prior_factors <- lapply(predictions, function(x) factor_gp_covariance(x$K, control$gp_jitter))
  state <- list(
    mean = do.call(cbind, lapply(predictions, `[[`, "mean")),
    mu = matrix(rep(vapply(fits, `[[`, numeric(1), "mean"), each = length(tt)), length(tt)),
    L = lapply(prior_factors, function(factor) t(factor$R)),
    bounds = attr(control, "state_bounds")
  )
  state$mean <- clip_state_to_bounds(state$mean, state$bounds)
  state$coordinates <- apply_state_bounds_to_coordinates(
    list(mu = state$mu, L = state$L, order = seq_len(ncol(state$mu))), state$bounds
  )
  state
}
