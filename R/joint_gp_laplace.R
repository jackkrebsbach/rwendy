# Experimental deterministic inference, separate from the fixed-lambda solver.
# theta = c(p, prior-whitened grid state); no Monte Carlo is used here.

.jgp_laplace_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop(name, " must be finite and positive.", call. = FALSE)
  x
}

.jgp_laplace_problem <- function(fit) {
  if (!inherits(fit, "jointgp") || is.null(fit$problem)) stop("Supply a solveWendyGP fit.")
  if (identical(fit$formulation,"latent"))
    stop("Laplace lambda inference currently requires the observed formulation.")
  if (!is.null(fit$problem$state_bounds) && any(fit$problem$state_bounds$bounded))
    stop("Laplace lambda inference does not yet integrate truncated state priors; state bounds are unsupported.")
  if (!isTRUE(fit$include_gp_prior))
    stop("Laplace lambda inference requires the GP prior.")
  if (!is.null(fit$beta) && !identical(as.numeric(fit$beta), 1))
    stop("Refit legacy tempered results before Laplace lambda inference.")
  if (!identical(fit$ode_weighting, "test_gram"))
    stop("Laplace lambda inference currently requires ode_weighting = 'test_gram'.")
  pr <- fit$problem; J <- pr$model$J; D <- pr$model$D; m <- nrow(pr$state$mu)
  k <- J + m * D
  lower <- c(pr$lower, rep(-Inf, m * D)); upper <- c(pr$upper, rep(Inf, m * D))
  if (length(lower) != k || length(upper) != k || anyNA(c(lower, upper)) || any(lower >= upper))
    stop("Invalid or zero-width parameter bounds.")
  W <- pr$weights$W; q <- nrow(W)
  if (is.null(q) || q < 1L || any(!is.finite(W))) stop("No finite retained weak space.")
  if (q < J) stop("The retained weak space cannot identify all flat-prior parameters.")
  # Complete the square in data likelihood + GP prior, without counting Y twice.
  R0 <- transform <- matrix(0, k, k)
  # Parameter columns/rows are exactly zero: there is NO parameter shrinkage.
  transform[cbind(seq_len(J), seq_len(J))] <- 1
  center <- c(fit$initial$p, numeric(m * D)); offset <- c(numeric(J), as.vector(pr$state$mu))
  for (d in seq_len(D)) {
    ii <- J + (d - 1L) * m + seq_len(m)
    A <- pr$H %*% pr$state$L[[d]] / fit$noise_sd[d]
    y <- (fit$Y[, d] - as.vector(pr$H %*% pr$state$mu[, d])) / fit$noise_sd[d]
    R <- chol(diag(m) + crossprod(A))
    R0[ii, ii] <- R; center[ii] <- .jgp_solve(R, crossprod(A, y))
    transform[ii, ii] <- pr$state$L[[d]]
  }
  evaluate <- function(theta, lambda, jacobian = TRUE) {
    point <- .jgp_unpack(theta, pr$state, J)
    v <- .jgp_weak_eval(point$U, point$p, pr$weak, pr$model, jacobian)
    w <- as.vector(W %*% v$r)
    jac <- NULL
    if (jacobian) {
      B <- cbind(v$Jp, v$Ju)
      for (d in seq_len(D)) {
        ii <- J + (d - 1L) * m + seq_len(m)
        B[, ii] <- B[, ii, drop = FALSE] %*% pr$state$L[[d]]
      }
      jac <- rbind(R0, sqrt(lambda) * W %*% B)
    }
    list(r = c(as.vector(R0 %*% (theta - center)), sqrt(lambda) * w), J = jac,
         S = sum(w^2))
  }
  list(evaluate = evaluate, center = center, start = pr$theta, lower = lower, upper = upper,
       q = q, transform = transform, offset = offset, J = J, D = D, m = m,
       parameter_prior = list(type = "flat", lower = pr$lower, upper = pr$upper,
         proper_on_finite_box = all(is.finite(c(pr$lower, pr$upper)))))
}

.jgp_laplace_control <- function(control = list()) {
  defaults <- list(maxit = 200L, refinements = 4L, hessian_step = 1e-4,
    quadrature_tol = .005, curvature_tol = .001, boundary_tol = .001,
    gradient_tol = 1e-5, tail_log_drop = 14, verbose = interactive())
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyDuplicated(names(control)) || any(!names(control) %in% names(defaults)))))
    stop("Invalid/unknown Laplace control.")
  defaults[names(control)] <- control; z <- defaults
  for (nm in setdiff(names(z), c("refinements", "verbose"))) .jgp_laplace_positive(z[[nm]], nm)
  for (nm in c("maxit", "refinements"))
    if (!is.numeric(z[[nm]]) || length(z[[nm]]) != 1L || !is.finite(z[[nm]]) ||
        z[[nm]] < 0 || z[[nm]] != floor(z[[nm]])) stop("Invalid ", nm, ".")
  if (!is.logical(z$verbose) || length(z$verbose) != 1L || is.na(z$verbose)) stop("Invalid verbose.")
  z
}

# Differentiate the full gradient, not just the residual Jacobian. These extra
# second-derivative terms affect both the marginal evidence and state covariance.
.jgp_laplace_hessian <- function(theta, evaluate, step) {
  k <- length(theta); H <- matrix(0, k, k)
  for (j in seq_len(k)) {
    h <- step * max(1, abs(theta[j])); plus <- minus <- theta
    plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
    vp <- evaluate(plus, TRUE); vm <- evaluate(minus, TRUE)
    H[, j] <- as.vector(crossprod(vp$J, vp$r) - crossprod(vm$J, vm$r)) / (2 * h)
  }
  (H + t(H)) / 2
}

.jgp_laplace_node <- function(problem, lambda, starts, control) {
  ev <- function(x, jacobian = TRUE) problem$evaluate(x, lambda, jacobian)
  ctl <- list(maxit = control$maxit, damping = 1e-3, gtol = control$gradient_tol / 10,
              ftol = 1e-13, xtol = 1e-11)
  candidates <- lapply(starts, function(x) tryCatch(
    .jgp_lm(x, ev, ctl, problem$lower, problem$upper), error = function(e) e))
  valid <- vapply(candidates, function(x) !inherits(x, "error") && is.finite(x$objective), logical(1))
  if (!any(valid)) stop("All conditional optimizations failed at lambda = ", lambda, ".")
  values <- vapply(candidates[valid], `[[`, numeric(1), "objective")
  best <- candidates[valid][[which.min(values)]]; theta <- best$par
  at <- ev(theta, TRUE)
  H1 <- .jgp_laplace_hessian(theta, ev, control$hessian_step)
  H <- .jgp_laplace_hessian(theta, ev, control$hessian_step / 2)
  R <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(R)) stop("Conditional full Hessian is not positive definite at lambda = ",
    signif(lambda, 6), "; Laplace integration stopped. No curvature ridge or dropped node was substituted.")
  C <- chol2inv(R)
  Q <- backsolve(R, diag(length(theta)))
  curvature_error <- norm(crossprod(Q, (H - H1) %*% Q), "F")
  # A union bound on Gaussian mass outside the parameter box. No truncated
  # Gaussian correction is silently invented for a mode close to a bound.
  s <- sqrt(pmax(0, diag(C)))
  outside <- sum(stats::pnorm((problem$lower - theta) / s) +
                   stats::pnorm((problem$upper - theta) / s, lower.tail = FALSE))
  gradient <- max(abs(as.vector(crossprod(at$J, at$r))) / sqrt(diag(H)))
  gap <- if (length(values) > 1L) max(values) - min(values) else NA_real_
  list(theta = theta, covariance = C, value = sum(at$r^2) / 2, S = at$S,
       logdet = 2 * sum(log(diag(R))), gradient = gradient, curvature_error = curvature_error,
       boundary_mass_bound = min(1, outside), start_objective_gap = gap,
       failed_starts = sum(!valid), gn_logdet = as.numeric(determinant(crossprod(at$J), logarithm = TRUE)$modulus))
}

.jgp_laplace_weights <- function(t, log_density) {
  dt <- diff(t); w <- c(dt[1L], head(dt, -1L) + tail(dt, -1L), tail(dt, 1L)) / 2
  scaled <- w * exp(log_density - max(log_density))
  scaled / sum(scaled)
}

.jgp_laplace_moments <- function(nodes, w, problem) {
  means <- vapply(nodes, function(x) as.vector(problem$transform %*% x$theta + problem$offset),
                  numeric(length(problem$start)))
  mean <- as.vector(means %*% w)
  conditional <- Reduce(`+`, Map(function(x, a) a * x$covariance, nodes, w))
  centered <- sweep(means, 1, mean, "-")
  covariance <- problem$transform %*% conditional %*% t(problem$transform) +
    tcrossprod(sweep(centered, 2, sqrt(w), "*"))
  list(mean = mean, covariance = covariance, sd = sqrt(pmax(0, diag(covariance))))
}

.jgp_laplace_integrate <- function(problem, lambda_prior, lambda_grid, control) {
  if (!is.list(lambda_prior) || length(lambda_prior) != 2L ||
      !setequal(names(lambda_prior), c("shape", "rate")))
    stop("Supply lambda_prior = list(shape = a, rate = b), a and b strictly positive.")
  a <- .jgp_laplace_positive(lambda_prior$shape, "lambda_prior$shape")
  b <- .jgp_laplace_positive(lambda_prior$rate, "lambda_prior$rate")
  if (is.null(lambda_grid)) {
    limits <- c(stats::qgamma(1e-7, shape = a, rate = b),
                stats::qgamma(1 - 1e-7, shape = a + problem$q / 2, rate = b))
    if (any(!is.finite(limits)) || any(limits <= 0)) stop("Supply a finite positive lambda_grid.")
    lambda_grid <- exp(seq(log(limits[1L]), log(limits[2L]), length.out = 25L))
  }
  if (!is.numeric(lambda_grid) || length(lambda_grid) < 5L || any(!is.finite(lambda_grid)) ||
      any(lambda_grid <= 0) || any(diff(lambda_grid) <= 0))
    stop("lambda_grid must contain at least five strictly increasing finite positive values.")
  cache <- new.env(parent = emptyenv()); solved <- numeric(); solutions <- list()
  node <- function(t) {
    key <- sprintf("%.17g", t)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
    starts <- list(problem$start, pmax(problem$lower, pmin(problem$upper, problem$center)))
    if (length(solved)) starts <- c(list(solutions[[which.min(abs(solved - t))]]), starts)
    ans <- .jgp_laplace_node(problem, exp(t), starts, control)
    # log density in t = log(lambda); the +t Jacobian is included.
    ans$log_density <- (a + problem$q / 2) * t - b * exp(t) - ans$value - ans$logdet / 2
    assign(key, ans, cache); solved <<- c(solved, t); solutions[[length(solutions) + 1L]] <<- ans$theta
    if (control$verbose) message("Laplace node lambda = ", signif(exp(t), 5), " complete.")
    ans
  }
  t <- log(lambda_grid)
  # Arithmetic construction of a grid plus explicit reference values can give
  # two effectively identical nodes (e.g. 100 and 100 - 1 ulp).
  keep <- c(TRUE, diff(t) > 64 * .Machine$double.eps *
              pmax(1, abs(head(t, -1L)), abs(tail(t, -1L))))
  t <- t[keep]
  if (length(t) < 5L) stop("Need five numerically distinct log-lambda nodes.")
  previous <- NULL; integration_error <- NA_real_
  history <- list()
  for (level in 0:control$refinements) {
    if (level > 0L) {
      mid <- (head(t, -1L) + tail(t, -1L)) / 2
      # First refinement spans the whole range. Later ones concentrate on
      # resolved posterior support; tail checks still apply to the full range.
      active <- if (level == 1L) rep(TRUE, length(mid)) else
        pmax(head(lp, -1L), tail(lp, -1L)) > max(lp) - control$tail_log_drop - 4
      t <- sort(c(t, mid[active]))
    }
    nodes <- lapply(t, node); lp <- vapply(nodes, `[[`, numeric(1), "log_density")
    w <- .jgp_laplace_weights(t, lp); moments <- .jgp_laplace_moments(nodes, w, problem)
    lambda_mean <- sum(w * exp(t)); lambda_sd <- sqrt(sum(w * (exp(t) - lambda_mean)^2))
    density <- exp(lp - max(lp))
    area <- diff(t) * (head(density, -1L) + tail(density, -1L)) / 2
    cdf <- c(0, cumsum(area) / sum(area))
    if (!is.null(previous)) integration_error <- max(
      abs(lambda_mean - previous$lambda_mean) / max(lambda_mean, .Machine$double.eps),
      abs(lambda_sd - previous$lambda_sd) / max(lambda_sd, .Machine$double.eps),
      abs(moments$mean - previous$moments$mean) / pmax(moments$sd, 1e-12),
      abs(moments$sd - previous$moments$sd) / pmax(moments$sd, 1e-12),
      abs(stats::approx(t, cdf, previous$t, rule = 2)$y - previous$cdf))
    previous <- list(lambda_mean = lambda_mean, lambda_sd = lambda_sd, moments = moments, t = t, cdf = cdf)
    history[[level + 1L]] <- data.frame(level = level, nodes = length(t),
      lambda_mean = lambda_mean, lambda_sd = lambda_sd, error = integration_error)
    if (is.finite(integration_error) && integration_error <= control$quadrature_tol) break
  }
  fields <- c("value", "S", "logdet", "gradient", "curvature_error", "boundary_mass_bound",
              "start_objective_gap", "failed_starts", "gn_logdet")
  table <- data.frame(lambda = exp(t), log_lambda = t, weight = w, log_density = lp)
  for (nm in fields) table[[nm]] <- vapply(nodes, `[[`, numeric(1), nm)
  # Endpoint density is a diagnostic, not a rigorous omitted-tail probability.
  tail_drop <- max(lp) - lp[c(1L, length(lp))]
  # A proper density need not have finite second moments under a flat p prior.
  # Check moment integrands, not only endpoint probability density.
  second <- vapply(nodes, function(x) pmax(diag(x$covariance) + x$theta^2,
    .Machine$double.xmin), numeric(length(problem$start)))
  moment_profiles <- rbind(sweep(log(second), 2, lp, "+"), lp + t, lp + 2 * t)
  moment_drop <- apply(moment_profiles, 1, max) - moment_profiles[, c(1L, length(lp)), drop = FALSE]
  moment_tail_drop <- apply(moment_drop, 2, min)
  bad <- table$gradient > control$gradient_tol | table$curvature_error > control$curvature_tol |
    table$boundary_mass_bound > control$boundary_tol | table$failed_starts > 0 |
    is.na(table$start_objective_gap) | table$start_objective_gap > .01
  unresolved_mass <- sum(w[bad])
  diagnostics <- list(quadrature_error = integration_error, endpoint_log_drop = tail_drop,
    moment_endpoint_log_drop = moment_tail_drop,
    flagged_node_mass = unresolved_mass, flagged_nodes = which(bad),
    numerical_passed = is.finite(integration_error) && integration_error <= control$quadrature_tol &&
      all(tail_drop >= control$tail_log_drop) && all(moment_tail_drop >= control$tail_log_drop) &&
      unresolved_mass <= control$boundary_tol,
    history = do.call(rbind, history), multimodality_ruled_out = FALSE,
    weak_quadrature_audited = FALSE)
  # Integrate a linearly interpolated log-lambda density for scalar quantiles;
  # moment calculations above use the same composite-trapezoid quadrature.
  quantile_t <- function(prob) {
    i <- min(length(t) - 1L, max(1L, findInterval(prob, cdf)))
    desired <- (prob - cdf[i]) * sum(area); h <- t[i + 1L] - t[i]
    if (desired <= 0) return(exp(t[i]))
    slope <- (density[i + 1L] - density[i]) / h
    advance <- if (abs(slope) * h < 1e-10 * max(density[i], 1e-300)) desired / density[i] else
      2 * desired / (density[i] + sqrt(max(0, density[i]^2 + 2 * slope * desired)))
    exp(t[i] + min(h, max(0, advance)))
  }
  list(lambda = c(mean = lambda_mean, sd = lambda_sd,
    stats::setNames(vapply(c(.025, .5, .975), quantile_t, numeric(1)), c("q025", "median", "q975"))),
    moments = moments, nodes = table, conditional = nodes, diagnostics = diagnostics,
    lambda_prior = list(shape = a, rate = b), retained_rank = problem$q)
}

#' Deterministic approximate joint inference of ODE precision, parameters and state
#'
#' @param fit A solveWendyGP fit with include_gp_prior = TRUE and
#'   ode_weighting = "test_gram". Its weak space and all GP/noise estimates stay fixed.
#' @param lambda_prior Proper Gamma prior list(shape = a, rate = b), with a,b > 0.
#' @param lambda_grid Increasing positive integration nodes. NULL initializes a
#'   log grid from prior quantiles; endpoint diagnostics may require a wider grid.
#' @param control List: maxit (200), refinements (4), hessian_step (1e-4),
#'   quadrature_tol (.005), curvature_tol (.001), boundary_tol (.001),
#'   gradient_tol (1e-5), tail_log_drop (14), verbose (interactive()). These are
#'   numerical controls, not model discrepancy or lambda calibration constants.
#' @details Conditional optimization uses multiple deterministic starts. The full
#'   Hessian is finite-differenced at two resolutions; no Hessian ridge is added.
#'   Laplace integration over (p,U) includes -log(det(H))/2. Log-lambda quadrature
#'   includes its change-of-variable Jacobian and the lambda^(q/2) constraint
#'   normalization, where q is the retained interior-plus-conditional-BL rank.
#'   The model uses a Gaussian constraint-observation factor, not a separately
#'   normalized dynamic prior on U given p. Conditional normals are combined by
#'   total expectation/covariance, so uncertainty in lambda is not discarded.
#'
#'   No sampling or Monte Carlo is performed. Integration is exact in (p,U) only
#'   for an unbounded Gaussian conditional problem; nonlinear problems use a
#'   local Laplace approximation. Bounds, unresolved competing modes, numerical
#'   curvature, quadrature resolution and lambda tails are diagnostic checks.
#'   A materially truncated conditional Gaussian is flagged, not corrected.
#'   Passing numerical checks is not validation of Laplace accuracy or coverage.
#'   GP hyperparameters, estimated noise and test selection uncertainty are not
#'   integrated. Weak quadrature/omitted directions need separate audits.
#'
#'   Finite lambda describes a soft weak constraint. Exactly zero discrepancy
#'   is a separate hard-constraint limit. Learning lambda does not directly
#'   minimize state or parameter MSE. No scientifically objective prior default
#'   or preferred lambda of 100 is assumed.
#' @return A jointgp_laplace list containing lambda marginal summaries, phat and
#'   U_hat (approximate marginal means), their joint covariance, integration nodes,
#'   conditional modes/covariances in GP-whitened coordinates, and diagnostics.
#'   ODE parameter priors are flat, restricted only by bounds already in fit.
#'   With unbounded flat priors, global posterior propriety is not guaranteed by
#'   positive local curvature and must be established separately.
#' @references Martino and Riebler (2019), Integrated Nested Laplace Approximations,
#'   https://arxiv.org/abs/1907.01248. This function is not an R-INLA implementation.
#' @export
inferWendyGPLambda <- function(fit, lambda_prior, lambda_grid = NULL, control = list()) {
  problem <- .jgp_laplace_problem(fit)
  control <- .jgp_laplace_control(control)
  ans <- .jgp_laplace_integrate(problem, lambda_prior, lambda_grid, control)
  J <- problem$J; m <- problem$m; D <- problem$D
  ans$phat <- ans$moments$mean[seq_len(J)]
  ans$U_hat <- matrix(ans$moments$mean[-seq_len(J)], m, D)
  ans$U_obs_hat <- fit$problem$H %*% ans$U_hat
  ans$U_sd <- matrix(ans$moments$sd[-seq_len(J)], m, D)
  ans$parameter_prior <- problem$parameter_prior
  ans$diagnostics$flat_prior_propriety_established <- problem$parameter_prior$proper_on_finite_box
  ans$bounds <- list(lower = fit$problem$lower, upper = fit$problem$upper)
  ans$tt <- fit$tt; ans$ode_weighting <- "test_gram"
  ans$conditioned_on <- c("GP hyperparameters", "measurement noise", "weak space", "weak metric")
  ans$control <- control; ans$call <- match.call()
  if (!ans$diagnostics$numerical_passed) warning(
    "Laplace numerical diagnostics are unresolved: inspect tails, grid resolution, bounds and node diagnostics. ",
    "This is not a validated posterior summary.", call. = FALSE)
  structure(ans, class = "jointgp_laplace")
}
