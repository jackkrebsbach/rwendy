# GP covariance, fitting, prediction, and grid-state preparation.

# The peak rescaling exp(eta) cancels under L2 normalization but avoids tiny
# numbers. The derivatives come from the same phi and derivative cache used
# by the existing weak-form code. Evaluation outside support is exactly zero.
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

factor_gp_covariance <- function(covariance, jitter = 0) {
  covariance <- (covariance + t(covariance)) / 2
  added <- jitter * max(diag(covariance))
  factor <- chol(covariance + diag(added, nrow(covariance)))
  list(R = factor, added = added)
}

solve_cholesky_system <- function(R, b) backsolve(R, forwardsolve(t(R), b))

evaluate_gp_radius <- function(x, coefficients, bounds) {
  # Preserve the location-by-coefficient Jacobian shape even for one location.
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

build_convolution_quadrature <- function(bounds, per_radius) {
  # Integration domain contains every admissible support for x in [0,1].
  left <- -bounds[2]
  right <- 1 + bounds[2]
  intervals <- ceiling((right - left) / (bounds[1] / per_radius))
  nodes <- seq(left, right, length.out = intervals + 1L)
  weights <- rep((right - left) / intervals, length(nodes))
  weights[c(1L, length(weights))] <- weights[c(1L, length(weights))] / 2
  list(nodes = nodes, weights = weights)
}

build_bump_features <- function(x, coefficients, bounds, quad, eta, derivatives = FALSE) {
  radius <- evaluate_gp_radius(x, coefficients, bounds)
  scaled_distance <- sweep(outer(x, quad$nodes, "-"), 1, radius$r, "/")
  bump <- evaluate_bump(scaled_distance, eta = eta)
  features <- sweep(bump, 2, sqrt(quad$weights), "*") / sqrt(radius$r)
  norms <- sqrt(rowSums(features^2))
  if (any(!is.finite(norms) | norms <= 0)) stop("Unresolved convolution quadrature.")
  normalized_features <- features / norms
  feature_derivative <- NULL
  if (derivatives) {
    raw_derivative <- sweep(
      -scaled_distance * evaluate_bump(scaled_distance, 1L, eta) - 0.5 * bump,
      2, sqrt(quad$weights), "*"
    ) / sqrt(radius$r)
    feature_derivative <- (raw_derivative - normalized_features * rowSums(normalized_features * raw_derivative)) / norms
  }
  list(B = normalized_features, dB = feature_derivative, dr = radius$derivative)
}

evaluate_gp_correlation <- function(x, coefficients, bounds, quad, control, derivatives = FALSE) {
  if (control$kernel == "bump") {
    features <- build_bump_features(x, coefficients, bounds, quad, control$bump_eta, derivatives)
    R <- tcrossprod(features$B)
    dR <- if (derivatives) {
      lapply(seq_along(coefficients), function(j) {
        part <- tcrossprod(features$dB * features$dr[, j], features$B)
        part + t(part)
      })
    } else {
      NULL
    }
    return(list(R = R, dR = dR))
  }
  # Paciorek-Schervish determinant factor and averaged squared local scales.
  radius <- evaluate_gp_radius(x, coefficients, bounds)
  r <- radius$r
  r2 <- r^2
  mean_r2 <- outer(r2, r2, "+") / 2
  z <- sqrt(5) * abs(outer(x, x, "-")) / sqrt(mean_r2)
  prefactor <- sqrt(outer(r, r) / mean_r2)
  decay <- exp(-z)
  R <- prefactor * (1 + z + z^2 / 3) * decay
  # radius$derivative is d log(radius) / d coefficient. Differentiate both the
  # determinant prefactor and the Matern factor through the averaged r^2.
  # This form never divides by distance or z, so repeated times and diagonal
  # entries are well defined (the unit-diagonal correlation has zero gradient).
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
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps) * max(1, abs(center))) {
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
  quad <- build_convolution_quadrature(bounds, control$gp_quad_per_radius)
  best <- NULL
  for (refine in 0:control$gp_quad_max_refine) {
    evaluate <- function(par) {
      coef <- par[seq_len(n_coef)]
      variance_scale <- exp(par[n_coef + 1L])
      correlation <- evaluate_gp_correlation(x, coef, bounds, quad, control, TRUE)
      covariance <- if (known_noise) {
        variance_scale * correlation$R + diag(noise2, n_obs)
      } else {
        correlation$R + diag(variance_scale, n_obs)
      }
      # An optimizer trial can produce a numerically indefinite covariance.
      # Reject that trial so the marginal-likelihood search can continue.
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
    starts <- if (!is.null(best)) {
      list(best$par)
    } else {
      lapply(seq_len(control$gp_restarts), function(j) {
        fraction <- seq(0.3, 0.8, length.out = control$gp_restarts)[j]
        c(stats::qlogis(fraction), rep(0, n_coef - 1L), if (known_noise) 0 else log(0.03))
      })
    }
    fits <- lapply(starts, function(init) {
      # Marginal likelihood includes a log determinant, so use a scalar solver.
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
    finer <- build_convolution_quadrature(bounds, control$gp_quad_per_radius * 2^(refine + 1L))
    fine_correlation <- evaluate_gp_correlation(x, best$par[seq_len(n_coef)], bounds, finer, control)$R
    quad_error <- max(abs(value$R - fine_correlation))
    if (control$kernel != "bump" || quad_error <= control$gp_quad_tol) break
    if (refine == control$gp_quad_max_refine) stop("GP convolution quadrature did not converge.")
    quad <- finer
  }
  observation_covariance <- value$tau2 * scale^2 * value$R + diag(value$noise2 * scale^2, n_obs)
  observation_factor <- factor_gp_covariance(observation_covariance)$R
  structure(list(
    tt = tt, x = x, y = y, origin = origin, span = span,
    mean = center + scale * value$mean, tau2 = scale^2 * value$tau2,
    noise2 = scale^2 * value$noise2, radius_coef = best$par[seq_len(n_coef)],
    radius_bounds = bounds, quad = quad, control = control, Ryy = observation_factor,
    alpha = solve_cholesky_system(observation_factor, y - center - scale * value$mean),
    convergence = best$convergence, message = best$message, objective = best$objective,
    quadrature_error = quad_error, quadrature_refinements = refine
  ), class = "wendygp_gp")
}

predict_component_gp <- function(fit, tt) {
  # Evaluate training and prediction locations with the same feature quadrature.
  x <- (tt - fit$origin) / fit$span
  all_locations <- c(fit$x, x)
  n_obs <- length(fit$x)
  prediction_indices <- n_obs + seq_along(x)
  R <- evaluate_gp_correlation(all_locations, fit$radius_coef, fit$radius_bounds, fit$quad, fit$control)$R
  cross_covariance <- fit$tau2 * R[prediction_indices, seq_len(n_obs), drop = FALSE]
  prior_covariance <- fit$tau2 * R[prediction_indices, prediction_indices, drop = FALSE]
  conditional_factor <- forwardsolve(t(fit$Ryy), t(cross_covariance))
  posterior_covariance <- prior_covariance - crossprod(conditional_factor)
  list(
    mean = as.vector(fit$mean + cross_covariance %*% fit$alpha), K = prior_covariance,
    Sigma = (posterior_covariance + t(posterior_covariance)) / 2
  )
}

gp_cross_covariance <- function(fit, t1, t2) {
  x <- (t1 - fit$origin) / fit$span
  y <- (t2 - fit$origin) / fit$span
  if (fit$control$kernel == "bump") {
    x_features <- build_bump_features(x, fit$radius_coef, fit$radius_bounds, fit$quad, fit$control$bump_eta)$B
    y_features <- build_bump_features(y, fit$radius_coef, fit$radius_bounds, fit$quad, fit$control$bump_eta)$B
    return(fit$tau2 * tcrossprod(x_features, y_features))
  }
  rx <- evaluate_gp_radius(x, fit$radius_coef, fit$radius_bounds)$r
  ry <- evaluate_gp_radius(y, fit$radius_coef, fit$radius_bounds)$r
  mean_r2 <- outer(rx^2, ry^2, "+") / 2
  z <- sqrt(5) * abs(outer(x, y, "-")) / sqrt(mean_r2)
  fit$tau2 * sqrt(outer(rx, ry) / mean_r2) * (1 + z + z^2 / 3) * exp(-z)
}

prepare_gp_state <- function(fits, tt, control) {
  predictions <- lapply(fits, predict_component_gp, tt = tt)
  prior_factors <- lapply(predictions, function(x) factor_gp_covariance(x$K, control$gp_jitter))
  state <- list(
    mean = do.call(cbind, lapply(predictions, `[[`, "mean")),
    mu = matrix(rep(vapply(fits, `[[`, numeric(1), "mean"), each = length(tt)), length(tt)),
    L = lapply(prior_factors, function(factor) t(factor$R)),
    Sigma = Map(function(prediction, factor) prediction$Sigma + diag(factor$added, length(tt)), predictions, prior_factors),
    jitter = vapply(prior_factors, `[[`, numeric(1), "added"),
    include_gp_prior = !isFALSE(control$include_gp_prior)
  )
  use_prior <- state$include_gp_prior
  state$bounds <- attr(control, "state_bounds")
  state$mean <- clip_state_to_bounds(state$mean, state$bounds)
  state$coordinates <- apply_state_bounds_to_coordinates(list(
    mu = if (use_prior) state$mu else state$mu * 0,
    L = if (use_prior) state$L else rep(list(1), ncol(state$mu)),
    order = seq_len(ncol(state$mu))
  ), state$bounds)
  state$prior_whitened <- use_prior && (is.null(state$bounds) || !any(state$bounds$bounded))
  state
}
